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

#ifndef XDMF_WRITER_BACKENDS_POSIX_H
#define XDMF_WRITER_BACKENDS_POSIX_H

#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

#include "Base.h"

namespace xdmfwriter
{

namespace backends
{

class Posix : public Base
{
private:
	/** File handles */
	int* m_fh;

public:
	Posix(unsigned int topoTypeSize)
		: Base(topoTypeSize),
		  m_fh(0L)
	{
	}

	virtual ~Posix()
	{
		close();
	}

	/**
	 *
	 * @param prefix
	 * @param variables List of variables which are written
	 * @param totalSize The total number of cells and vertices
	 * @param localSize The local number of cells and vertices
	 * @param offset The offset of the local process for cells and vertices
	 * @param cells Cell data (only required for <code>create == true</code>)
	 * @param vertices Vertex data (only required for <code>create == true</code>)
	 * @param partition Partition data (only required for <code>create == true</code>)
	 * @param create
	 */
	void open(const std::string &prefix, const std::vector<const char*> &variables,
			const unsigned long totalSize[2], const unsigned int localSize[2],
			const unsigned long offset[2],
			const unsigned long* cells, const double* vertices,
			const unsigned int* partition, bool create = true)
	{
		Base::open(variables.size(), totalSize[0], localSize[0], offset[0]);

		if (create) {
			// Create backups
			backup(filename(prefix, "connect"));
			backup(filename(prefix, "geometry"));
			backup(filename(prefix, "partition"));
			for (std::vector<const char*>::const_iterator i = variables.begin();
					i != variables.end(); i++)
				backup(filename(prefix, *i));

#ifdef PARALLEL
			// Make sure the file is moved before continuing
			MPI_Barrier(comm());
#endif // PARALLEL

			// Create connect file
			int fh = open(filename(prefix, "connect").c_str());
			write(fh, cells, offset[0]*m_topoTypeSize*sizeof(unsigned long),
					localSize[0]*m_topoTypeSize*sizeof(unsigned long));
			checkErr(::close(fh));

			// Create geometry file
			fh = open(filename(prefix, "geometry").c_str());
			write(fh, vertices, offset[1]*3*sizeof(double),
					localSize[1]*3*sizeof(double));
			checkErr(::close(fh));

			// Create partition file
			fh = open(filename(prefix, "partition").c_str());
			write(fh, partition, offset[0]*sizeof(unsigned int),
					localSize[0]*sizeof(unsigned int));
			checkErr(::close(fh));
		}

		// Create/open variable files
		m_fh = new int[variables.size()];
		for (unsigned int i = 0; i < numVars(); i++)
			m_fh[i] = open(filename(prefix, variables[i]).c_str());
	}

	void writeData(unsigned int timestep, unsigned int id, const double* data)
	{
		size_t offset = (totalCells() * timestep + offsetCells()) * sizeof(double);
		write(m_fh[id], data, offset, localCells() * sizeof(double));
	}

	void flush()
	{
		for (unsigned int i = 0; i < numVars(); i++)
			checkErr(fsync(m_fh[i]));
	}

	/**
	 * Closes the HDF5 file (should be done before MPI_Finalize is called)
	 */
	void close()
	{
		if (m_fh) {
			for (unsigned int i = 0; i < numVars(); i++)
				checkErr(::close(m_fh[i]));

			delete [] m_fh;
			m_fh = 0L; // Indicates closed files
		}
	}

public:
	static const char* format()
	{
		return "Binary";
	}

	/**
	 * @param prefix
	 * @param var
	 * @return The location name in the XDMF file for a data item
	 */
	static std::string dataItemLocation(const char* prefix, const char* var)
	{
		return filename(prefix, var);
	}

private:
	/**
	 * @return The file name for a given variable
	 */
	static std::string filename(const std::string &prefix, const char* var)
	{
		return prefix + "_" + var + ".bin";
	}

	static int open(const char* filename)
	{
		int fh = open64(filename, O_WRONLY | O_CREAT,
				S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
		checkErr(fh);
		return fh;
	}

	static void write(int fh, const void* buffer, size_t offset, size_t size)
	{
		checkErr(lseek64(fh, offset, SEEK_SET));

		const char* buf = reinterpret_cast<const char*>(buffer);
		while (size > 0) {
			ssize_t written = ::write(fh, buf, size);
			if (written <= 0)
				checkErr(written, size);
			buf += written;
			size -= written;
		}
	}

	template<typename T>
	static void checkErr(T ret)
	{
		if (ret < 0)
			logError() << "An POSIX error occurred in the XDMF writer:" << strerror(errno);
	}

	/**
	 * Can be used to check read/write errors
	 *
	 * @param ret The return value
	 * @param target The expected return value (> 0)
	 */
	template<typename T, typename U>
	static void checkErr(T ret, U target)
	{
		checkErr(ret);
		if (ret != static_cast<T>(target))
			logError() << "Error in XDMF writer:"
				<< target << "bytes expected;" << ret << "bytes gotten";
	}
};

}

}

#endif // XDMF_WRITER_BACKENDS_POSIX_H
