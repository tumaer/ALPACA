/**
 * @file
 *  This file is part of XdmfWriter
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2014-2015, Technische Universitaet Muenchen.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from this
 *     software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef BLOCK_BUFFER_H
#define BLOCK_BUFFER_H

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <stdint.h>

#include "submodules/utils/logger.h"

/**
 * Splits a buffer into equal sized blocks and distributes the
 * blocks across all ranks
 */
class BlockBuffer
{
private:
	/** The communicator that should be used */
	MPI_Comm m_comm;

	/** Our rank */
	int m_rank;

	/** The basic data type */
	MPI_Datatype m_type;

	/** Size of the data type */
	MPI_Aint m_extent;

	/** Number of elements before the exchange */
	unsigned int m_inCount;

	/** Rank to which we send our data */
	int m_sendRank;

	/** Number of elements we send to the left */
	unsigned int m_sendCount;

	/** Number of ranks from which we receive */
	unsigned int m_recvRanks;

	/** Number of elements we receive from the right */
	unsigned int *m_recvCounts;

	/** Number of elements we will have after and exchange */
	unsigned int m_outCount;

	/** Requests used for sending/receiving (all allocate them once) */
	MPI_Request *m_requests;

public:
	BlockBuffer()
		: m_comm(MPI_COMM_NULL), m_rank(-1), m_type(MPI_DATATYPE_NULL), m_extent(0), m_inCount(0),
		  m_sendRank(0), m_sendCount(0), m_recvRanks(0), m_recvCounts(0L), m_outCount(1), m_requests(0L)
	{
	}

	virtual ~BlockBuffer()
	{
		delete [] m_recvCounts;
	}

	/**
	 * @param comm The MPI communicator
	 * @param inCount The number of elements in the local buffer
	 * @param dataType The MPI data type of the elements
	 * @param blockSize The size of the blocks that should be created (in bytes)
	 */
	void init(MPI_Comm comm, unsigned int inCount, MPI_Datatype dataType, unsigned long blockSize)
	{
		m_comm = comm;
		MPI_Comm_rank(comm, &m_rank);

		MPI_Aint lb;
		MPI_Type_get_extent(dataType, &lb, &m_extent);
		m_type = dataType;

		if (blockSize % m_extent != 0)
			logError() << "The block size must be a multiple of the data type size";

		unsigned long localSize = inCount * m_extent;
		unsigned long localStart;
		MPI_Scan(&localSize, &localStart, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);
		localStart -= localSize;

		// Determine what we need to send, what is ours and what we receive
		unsigned long sendSize = (blockSize - localStart % blockSize) % blockSize;
		sendSize = std::min(sendSize, localSize);
		unsigned long recvSize;
		if (sendSize >= localSize)
			recvSize = 0;
		else
			recvSize = (blockSize - (localStart + localSize) % blockSize) % blockSize;
		assert((localSize - sendSize + recvSize) % blockSize == 0);
		unsigned int localBlocks = (localSize - sendSize + recvSize) / blockSize;

		m_sendCount = sendSize / m_extent;

		// Compute all the ranks, from which we receive and to which we send
		int numProcs;
		MPI_Comm_size(comm, &numProcs);

		unsigned int *blocks = new unsigned int[numProcs];
		MPI_Allgather(&localBlocks, 1, MPI_UNSIGNED, blocks, 1, MPI_UNSIGNED, comm);

		m_sendRank = m_rank-1;
		if (m_rank > 0) {
			// Only if we have a left process
			while (blocks[m_sendRank] == 0) {
				m_sendRank--;
				assert(m_sendRank >= 0);
			}
		}

		int lastRecvRank = m_rank+1;
		if (localBlocks == 0 || m_rank == numProcs-1)
			lastRecvRank = m_rank;
		else {
			// Only if we receive any elements
			while (blocks[lastRecvRank] == 0 && lastRecvRank+1 < numProcs)
				lastRecvRank++;
		}
		m_recvRanks = lastRecvRank - m_rank;

		delete [] blocks;

		// Exchange the elements we transfer from/to each rank
		m_requests = new MPI_Request[m_recvRanks + 1];
		m_recvCounts = new unsigned int[m_recvRanks];
		if (m_rank == 0)
			m_requests[0] = MPI_REQUEST_NULL;
		else
			MPI_Isend(&m_sendCount, 1, MPI_UNSIGNED, m_sendRank, 0,
					comm, &m_requests[0]);
		for (unsigned int i = 0; i < m_recvRanks; i++)
			MPI_Irecv(&m_recvCounts[i], 1, MPI_UNSIGNED, m_rank+i+1, 0,
					comm, &m_requests[i+1]);
		MPI_Waitall(m_recvRanks+1, m_requests, MPI_STATUSES_IGNORE); // Important: This sets all requests to MPI_REQEUST_NULL

		// Compute how many elements we are responsible for
		m_inCount = inCount;
		m_outCount = inCount - m_sendCount;
		for (unsigned int i = 0; i < m_recvRanks; i++)
			m_outCount += m_recvCounts[i];
	}

	/**
	 * @return <code>False</code> before initialization was called, <code>true</code> afterwards
	 */
	bool isInitialized() const
	{
		return m_rank >= 0;
	}

	/**
	 * Before initialization this will always return 1.
	 *
	 * @return The number of elements after the exchange
	 */
	unsigned int count() const
	{
		return m_outCount;
	}

	/**
	 * Exchanges the data according to the configuration
	 *
	 * @in The local data before the exchange.
	 * @out The local data after the exchange. The caller is responsible
	 *  for allocating this buffer.
	 */
	void exchange(const void* in, void* out)
	{
		_exchange(in, m_type, m_extent, 1, out);
	}

	/**
	 * Does the same as <code>exchange(void*, void*)</code> but creates block of
	 * a continues type <code>type</code> with <code>count</code> elements.
	 * The resulting block have the size type * count / orig_type * orig_block_size.
	 */
	void exchange(const void* in, MPI_Datatype type, unsigned int count, void* out)
	{
		MPI_Aint lb;
		MPI_Aint extent;
		MPI_Type_get_extent(type, &lb, &extent);

		_exchange(in, type, extent, count, out);
	}

	/**
	 * Exchanges any number of elements but does not garantie a specific block size. Makes
	 * sure that only processes with <code>count() > 0</code> get elements.
	 *
	 * @param count Number of variables of type <code>type</code> for one element
	 * @param elemIn Number of elements on this processes
	 * @param[out] elemOut Number of elements this processes has after the exchange
	 * @return Buffer the elements after the exchange. Use <code>free()</code>
	 *  to free the memory
	 */
	template<typename T>
	T* exchangeAny(const T* in, MPI_Datatype type, unsigned int count,
			unsigned int elemIn, unsigned int &elemOut)
	{
		MPI_Aint lb;
		MPI_Aint extent;
		MPI_Type_get_extent(type, &lb, &extent);

		extent *= count;

		// Recv size
		unsigned int *recvCounts = new unsigned int[m_recvRanks];
		memset(recvCounts, 0, m_recvRanks*sizeof(unsigned int));
		for (unsigned int i = 0; i < m_recvRanks; i++) {
			if (m_recvCounts[i] > 0)
				MPI_Irecv(&recvCounts[i], 1, MPI_UNSIGNED, m_rank+i+1, 0,
						m_comm, &m_requests[i+1]);
		}

		// Send size
		unsigned int sendCount = (m_outCount > 0 ? 0 : elemIn); // Send elements if we don't own any
		if (m_sendCount > 0)
			MPI_Isend(&sendCount, 1, MPI_UNSIGNED, m_sendRank, 0, m_comm, &m_requests[0]);

		MPI_Waitall(m_recvRanks+1, m_requests, MPI_STATUSES_IGNORE);

		// How many elements are we responsible for?
		elemOut = elemIn - sendCount;
		for (unsigned int i = 0; i < m_recvRanks; i++)
			elemOut += recvCounts[i];

		uint8_t *buf = 0L;
		if (elemOut > 0)
			buf = new uint8_t[elemOut * extent];

		// Recv elements
		unsigned int elemOffset = elemIn;
		for (unsigned int i = 0; i < m_recvRanks; i++) {
			if (recvCounts[i] > 0) {
				MPI_Irecv(buf+elemOffset*extent, recvCounts[i]*count, type,
						m_rank+i+1, 0, m_comm, &m_requests[i+1]);
				elemOffset += recvCounts[i];
			}
		}

		// Send elements or copy
		if (sendCount > 0)
			MPI_Isend(const_cast<T*>(in), sendCount*count, type, m_sendRank, 0,
					m_comm, &m_requests[0]);
		else
			memcpy(buf, in, elemIn*extent);

		MPI_Waitall(m_recvRanks+1, m_requests, MPI_STATUSES_IGNORE);

		return reinterpret_cast<T*>(buf);
	}

private:
	/**
	 * @param extent Extent of a single element
	 * @param count Number of elements
	 */
	void _exchange(const void* in, MPI_Datatype type, MPI_Aint extent, unsigned int count, void* out)
	{
		extent *= count;

		// Recv elements
		unsigned int elemOffset = m_inCount-m_sendCount;
		for (unsigned int i = 0; i < m_recvRanks; i++) {
			if (m_recvCounts[i] > 0) {
				MPI_Irecv(static_cast<uint8_t*>(out)+elemOffset*extent,
						m_recvCounts[i]*count, type, m_rank+i+1, 0, m_comm, &m_requests[i+1]);
				elemOffset += m_recvCounts[i];
			}
		}

		// Send first elements
		if (m_sendCount > 0)
			MPI_Isend(const_cast<void*>(in), m_sendCount*count, type, m_sendRank, 0,
					m_comm, &m_requests[0]);

		// Copy local elements
		memcpy(out, static_cast<const uint8_t*>(in)+m_sendCount*extent,
				(m_inCount-m_sendCount)*extent);

		MPI_Waitall(m_recvRanks+1, m_requests, MPI_STATUSES_IGNORE);
	}

public:
	static void free(void *buf)
	{
		delete [] static_cast<uint8_t*>(buf);
	}

};

#endif // BLOCK_BUFFER_H
