/**
 * @file
 *  This file is part of XdmfWriter
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2016, Technische Universitaet Muenchen.
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

#ifndef XDMF_WRITER_BACKENDS_BASE_H
#define XDMF_WRITER_BACKENDS_BASE_H

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

#include <string>

#include "3rdParty/xdmf_writer/submodules/utils/logger.h"

namespace xdmfwriter
{

namespace backends
{

class Base
{
protected:
	const unsigned int m_topoTypeSize;

private:
#ifdef PARALLEL
	MPI_Comm m_comm;
#endif // PARALLEL

	/** Rank of this process */
	int m_rank;

	/** The number of variables */
	unsigned int m_numVars;

	/** Total number of cells */
	unsigned long m_totalCells;
	/** Local number of cells */
	unsigned int m_localCells;
	/** Offset where we start writing our part */
	unsigned long m_offsetCells;

protected:
	Base(unsigned int topoTypeSize)
		: m_topoTypeSize(topoTypeSize),
#ifdef PARALLEL
		  m_comm(MPI_COMM_WORLD),
#endif // PARALLEL
		  m_rank(0),
		  m_numVars(0),
		  m_totalCells(0),
		  m_localCells(0),
		  m_offsetCells(0)
	{
	}

public:
#ifdef PARALLEL
	/**
	 * Sets the communicator that should be used. Default is MPI_COMM_WORLD.
	 */
	void setComm(MPI_Comm comm)
	{
		m_comm = comm;
		MPI_Comm_rank(comm, &m_rank);
	}
#endif // PARALLEL

protected:
	void open(unsigned int numVars, unsigned long totalCells,
			unsigned int localCells, unsigned long offsetCells)
	{
		m_numVars = numVars;
		m_totalCells = totalCells;
		m_localCells = localCells;
		m_offsetCells = offsetCells;
	}

	/**
	 * Backup an existing backend file
	 */
	void backup(const std::string &file)
	{
		// Backup any existing file
		struct stat statBuffer;
		if (m_rank == 0 && stat(file.c_str(), &statBuffer) == 0) {
			logWarning() << "File" << file << "already exists. Creating backup.";
			rename(file.c_str(), (file + ".bak").c_str());
		}
	}

#ifdef PARALLEL
	MPI_Comm comm() const
	{
		return m_comm;
	}
#endif // PARALLEL

	unsigned int numVars() const
	{
		return m_numVars;
	}

	unsigned long totalCells() const
	{
		return m_totalCells;
	}

	unsigned int localCells() const
	{
		return m_localCells;
	}

	unsigned long offsetCells() const
	{
		return m_offsetCells;
	}
};

}

}

#endif // XDMF_WRITER_BACKENDS_BASE_H
