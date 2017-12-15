/**
 * @file
 *  This file is part of XdmfWriter
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2014-2016, Technische Universitaet Muenchen.
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

#ifndef PARALLEL_VERTEX_FILTER_H
#define PARALLEL_VERTEX_FILTER_H

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <stdint.h> // TODO switch to cstdint as soon as all builds enable C++11
#include <vector>

#include "submodules/utils/logger.h"

#include "scorep_wrapper.h"

/**
 * Filters duplicate vertices in parallel
 */
class ParallelVertexFilter
{
private:
	/**
	 * Compares 3D-vertex indices according to the vertices
	 */
	class IndexedVertexComparator
	{
	private:
		const double *m_vertices;

	public:
		IndexedVertexComparator(const double *vertices)
			: m_vertices(vertices)
		{
		}

		bool operator() (unsigned int i, unsigned int j)
		{
			i *= 3;
			j *= 3;

			return (m_vertices[i] < m_vertices[j])
					|| (m_vertices[i] == m_vertices[j] && m_vertices[i+1] < m_vertices[j+1])
					|| (m_vertices[i] == m_vertices[j] && m_vertices[i+1] == m_vertices[j+1]
						&& m_vertices[i+2] < m_vertices[j+2]);
		}
	};

private:
	/** The communicator we use */
	MPI_Comm m_comm;

	/** Our rank */
	int m_rank;

	/** #Processes */
	int m_numProcs;

	/** Global id after filtering */
	unsigned long *m_globalIds;

	/** Number of local vertices after filtering */
	unsigned int m_numLocalVertices;

	/** Local vertices after filtering */
	double *m_localVertices;

	/** MPI data type consisting of three doubles */
	MPI_Datatype m_vertexType;

public:
	ParallelVertexFilter(MPI_Comm comm = MPI_COMM_WORLD)
		: m_comm(comm), m_globalIds(0L), m_numLocalVertices(0), m_localVertices(0L)
	{
		MPI_Comm_rank(comm, &m_rank);
		MPI_Comm_size(comm, &m_numProcs);

		MPI_Type_contiguous(3, MPI_DOUBLE, &m_vertexType);
		MPI_Type_commit(&m_vertexType);
	}

	virtual ~ParallelVertexFilter()
	{
		delete [] m_globalIds;
		delete [] m_localVertices;

		MPI_Type_free(&m_vertexType);
	}

	/**
	 * @param vertices Vertices that should be filtered, must have the size <code>numVertices * 3</code>
	 */
	void filter(unsigned int numVertices, const double *vertices)
	{
		SCOREP_USER_REGION("ParallelVertexFilter_Filter", SCOREP_USER_REGION_TYPE_FUNCTION);

		// Chop the last 4 bits to avoid numerical errors
		double *roundVertices = new double[numVertices*3];
		removeRoundError(vertices, numVertices*3, roundVertices);

		// Create indices and sort them locally
		unsigned int *sortIndices = new unsigned int[numVertices];
		createSortedIndices(roundVertices, numVertices, sortIndices);

		// Select BUCKETS_PER_RANK-1 splitter elements
		double localSplitters[BUCKETS_PER_RANK-1];
#if 0 // Use omp only if we create a larger amount of buckets
#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
#endif
		for (int i = 0; i < BUCKETS_PER_RANK-1; i++) {
			unsigned long vrtxIndex = static_cast<unsigned long>(i)
					* static_cast<unsigned long>(numVertices)
					/ static_cast<unsigned long>(BUCKETS_PER_RANK-1);
			assert(vrtxIndex < numVertices);

			localSplitters[i] = roundVertices[sortIndices[vrtxIndex]*3];
		}

		// Collect all splitter elements on rank 0
		double *allSplitters = 0L;

		if (m_rank == 0)
			allSplitters = new double[m_numProcs * (BUCKETS_PER_RANK-1)];

		MPI_Gather(localSplitters, BUCKETS_PER_RANK-1, MPI_DOUBLE,
				allSplitters, BUCKETS_PER_RANK-1, MPI_DOUBLE,
				0, m_comm);

		// Sort splitter elements
		if (m_rank == 0)
			std::sort(allSplitters, allSplitters + (m_numProcs * (BUCKETS_PER_RANK-1)));

		// Distribute splitter to all processes
		double *splitters = new double[m_numProcs-1];

		if (m_rank == 0) {
#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
#endif
			for (int i = 0; i < m_numProcs-1; i++) {
				unsigned long spltIndex = (i+1) * (BUCKETS_PER_RANK-1);
				assert(spltIndex < static_cast<unsigned int>(m_numProcs * (BUCKETS_PER_RANK-1)));

				splitters[i] = allSplitters[spltIndex];
			}
		}

		MPI_Bcast(splitters, m_numProcs-1, MPI_DOUBLE, 0, m_comm);

		delete [] allSplitters;

		// Determine the bucket for each vertex
		unsigned int *bucket = new unsigned int[numVertices];

#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
		for (unsigned int i = 0; i < numVertices; i++) {
			double* ub = std::upper_bound(splitters, splitters+m_numProcs-1, roundVertices[i*3]);

			bucket[i] = ub-splitters;
		}

		delete [] roundVertices;
		delete [] splitters;

		// Determine the (local and total) bucket size
		int *bucketSize = new int[m_numProcs];
		memset(bucketSize, 0, sizeof(int)*m_numProcs);
		for (unsigned int i = 0; i < numVertices; i++)
			bucketSize[bucket[i]]++;

		delete [] bucket;

		// Tell all processes what we are going to send them
		int *recvSize = new int[m_numProcs];

		MPI_Alltoall(bucketSize, 1, MPI_INT, recvSize, 1, MPI_INT, m_comm);

		unsigned int numSortVertices = 0;
#ifdef _OPENMP
		#pragma omp parallel for schedule(static) reduction(+: numSortVertices)
#endif
		for (int i = 0; i < m_numProcs; i++)
			numSortVertices += recvSize[i];

		// Create sorted send buffer
		double *sendVertices = new double[3 * numVertices];
#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
		for (unsigned int i = 0; i < numVertices; i++) {
			memcpy(&sendVertices[i*3], &vertices[sortIndices[i]*3], sizeof(double)*3);
		}

		// Allocate buffer for the vertices and exchange them
		double *sortVertices = new double[3 * numSortVertices];

		int *sDispls = new int[m_numProcs];
		int *rDispls = new int[m_numProcs];
		sDispls[0] = 0;
		rDispls[0] = 0;
		for (int i = 1; i < m_numProcs; i++) {
			sDispls[i] = sDispls[i-1] + bucketSize[i-1];
			rDispls[i] = rDispls[i-1] + recvSize[i-1];
		}
		MPI_Alltoallv(sendVertices, bucketSize, sDispls, m_vertexType, sortVertices, recvSize, rDispls, m_vertexType, m_comm);

		delete [] sendVertices;

		// Chop the last 4 bits to avoid numerical errors
		roundVertices = new double[numSortVertices*3];
		removeRoundError(sortVertices, numSortVertices*3, roundVertices);

		// Create indices and sort them (such that the vertices are sorted)
		unsigned int *sortSortIndices = new unsigned int[numSortVertices];
		createSortedIndices(roundVertices, numSortVertices, sortSortIndices);

		// Initialize the global ids we send back to the other processors
		unsigned long *gids = new unsigned long[numSortVertices];

		if (numSortVertices > 0) {
			gids[sortSortIndices[0]] = 0;
			for (unsigned int i = 1; i < numSortVertices; i++) {
				if (equals(&roundVertices[sortSortIndices[i-1]*3], &roundVertices[sortSortIndices[i]*3]))
					gids[sortSortIndices[i]] = gids[sortSortIndices[i-1]];
				else
					gids[sortSortIndices[i]] = gids[sortSortIndices[i-1]] + 1;
			}
		}

		delete [] roundVertices;

		// Create the local vertices list
		if (numSortVertices > 0)
			m_numLocalVertices = gids[sortSortIndices[numSortVertices-1]] + 1;
		else
			m_numLocalVertices = 0;
		delete [] m_localVertices;
		m_localVertices = new double[m_numLocalVertices * 3];
		for (unsigned int i = 0; i < numSortVertices; i++)
			memcpy(&m_localVertices[gids[i]*3], &sortVertices[i*3], sizeof(double)*3);

		delete [] sortVertices;

		// Get the vertices offset
		unsigned int offset = m_numLocalVertices;
		MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, m_comm);
		offset -= m_numLocalVertices;

		// Add offset to the global ids
#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
		for (unsigned int i = 0; i < numSortVertices; i++)
			gids[i] += offset;

		// Send result back
		unsigned long *globalIds = new unsigned long[numVertices];
		MPI_Alltoallv(gids, recvSize, rDispls, MPI_UNSIGNED_LONG,
				globalIds, bucketSize, sDispls, MPI_UNSIGNED_LONG, m_comm);

		delete [] bucketSize;
		delete [] recvSize;
		delete [] sDispls;
		delete [] rDispls;
		delete [] gids;

		// Assign the global ids to the correct vertices
		delete [] m_globalIds;
		m_globalIds = new unsigned long[numVertices];
#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
		for (unsigned int i = 0; i < numVertices; i++)
			m_globalIds[sortIndices[i]] = globalIds[i];

		delete [] sortIndices;
		delete [] globalIds;
	}

	/**
	 * @return The list of the global identifiers after filtering
	 */
	const unsigned long* globalIds() const
	{
		return m_globalIds;
	}

	/**
	 * @return Number of vertices this process is responsible for after filtering
	 */
	unsigned int numLocalVertices() const
	{
		return m_numLocalVertices;
	}

	/**
	 * @return The list of vertices this process is responsible for after filtering
	 */
	const double* localVertices() const
	{
		return m_localVertices;
	}

private:
	/**
	 * Removes round errors of double values by setting the last 4 bits
	 * (of the significand) to zero.
	 *
	 * @warning Only works if <code>value</code> ist not nan or infinity
	 * @todo This should work for arbitrary precision
	 */
	static double removeRoundError(double value)
	{
		static const uint64_t mask = ~0xFF;

		union FloatUnion {
			double f;
			uint64_t bits;
		};

		FloatUnion result;
		result.f = value;

		result.bits &= mask;

		return result.f;
	}

	/**
	 * Removes the round errors using {@link removeRoundError(double)}
	 *
	 * @param values The list of floating point values
	 * @param count Number of values
	 * @param[out] roundValues The list of rounded values
	 *  (the caller is responsible for allocating the memory)
	 */
	static void removeRoundError(const double *values, unsigned int count, double* roundValues)
	{
#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
		for (unsigned int i = 0; i < count; i++)
			roundValues[i] = removeRoundError(values[i]);
	}

	/**
	 * Creates the list of sorted indices for the vertices.
	 * The caller is responsible for allocating the memory.
	 */
	static void createSortedIndices(const double *vertices, unsigned int numVertices,
			unsigned int *sortedIndices)
	{

#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
		for (unsigned int i = 0; i < numVertices; i++)
			sortedIndices[i] = i;

		IndexedVertexComparator comparator(vertices);
		std::sort(sortedIndices, sortedIndices+numVertices, comparator);
	}

	/**
	 * Compares to vertices for equality
	 * Assumes that the rounding errors are removed.
	 */
	static bool equals(const double* vertexA, const double* vertexB)
	{
		return vertexA[0] == vertexB[0]
		       && vertexA[1] == vertexB[1]
		       && vertexA[2] == vertexB[2];
	}

	/** The total buckets we create is <code>BUCKETS_PER_RANK * numProcs</code> */
	const static int BUCKETS_PER_RANK = 8;
};

#endif // PARALLEL_VERTEX_FILTER_H
