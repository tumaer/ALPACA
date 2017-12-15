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

#ifndef XDMF_WRITER_BACKENDS_HDF5_H
#define XDMF_WRITER_BACKENDS_HDF5_H

#include <string>

#include <hdf5.h>
//#include "/global/HDF5/include/H5FDmpio.h"
#include <H5FDmpio.h>


#include "3rdParty/xdmf_writer/submodules/utils/env.h"

#include "Base.h"
#ifdef PARALLEL
#include "MPIInfo.h"
#endif // PARALLEL

namespace xdmfwriter
{

namespace backends
{

class HDF5 : public Base
{
private:
	hid_t m_hdfFile;

	hid_t m_hdfAccessList;

	/** Variable identifiers */
	hid_t *m_hdfVars;

	/** Memory space for writing one time step */
	hid_t m_hdfVarMemSpace;

public:
	HDF5(unsigned int topoTypeSize)
		: Base(topoTypeSize),
		  m_hdfFile(0),
          m_hdfAccessList(H5P_DEFAULT),
		  m_hdfVars(0L),
		  m_hdfVarMemSpace(0)
	{
	}

	virtual ~HDF5()
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

		hid_t h5plist = H5Pcreate(H5P_FILE_ACCESS);
		checkH5Err(h5plist);
		checkH5Err(H5Pset_libver_bounds(h5plist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
		checkH5Err(H5Pset_meta_block_size(h5plist, 1024*1024));
		hsize_t align = utils::Env::get<hsize_t>("XDMFWRITER_ALIGNMENT", 0);
		if (align > 0)
			checkH5Err(H5Pset_alignment(h5plist, 1, align));
#ifdef PARALLEL
		checkH5Err(H5Pset_fapl_mpio(h5plist, comm(), MPIInfo::get()));
#endif // PARALLEL

		// Assemble filename
		std::string filename = prefix + fileExention();

		// Create/open the file
		if (create) {
			// Backup existing file
			backup(filename);
#ifdef PARALLEL
			// Make sure the file is moved before continuing
			MPI_Barrier(comm());
#endif // PARALLEL

            m_hdfFile = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5plist);
		} else
			m_hdfFile = H5Fopen(filename.c_str(), H5F_ACC_RDWR, h5plist);
		checkH5Err(m_hdfFile);

		checkH5Err(H5Pclose(h5plist));

		// Parallel access for all datasets
#ifdef PARALLEL
		m_hdfAccessList = H5Pcreate(H5P_DATASET_XFER);
		checkH5Err(m_hdfAccessList);
		checkH5Err(H5Pset_dxpl_mpio(m_hdfAccessList, H5FD_MPIO_COLLECTIVE));
#endif // PARALLEL

		// Allocate array for the variable ids
		m_hdfVars = new hid_t[numVars()];

		if (create) {
			// Create connect dataset
			hsize_t connectDims[2] = {totalSize[0], m_topoTypeSize};
			hid_t h5connectSpace = H5Screate_simple(2, connectDims, 0L);
			checkH5Err(h5connectSpace);
			hid_t h5connect = H5Dcreate(m_hdfFile, "/connect", H5T_STD_U64LE, h5connectSpace,
					H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			checkH5Err(h5connect);

			// Create geometry dataset
			hsize_t geometryDims[2] = {totalSize[1], 3};
			hid_t h5geometrySpace = H5Screate_simple(2, geometryDims, 0L);
			checkH5Err(h5geometrySpace);
			hid_t h5geometry = H5Dcreate(m_hdfFile, "/geometry", H5T_IEEE_F64LE, h5geometrySpace,
					H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			checkH5Err(h5geometry);

                // NH 2016-10-13 - I HAVE NO IDEA WHAT I'M DOING
//                hsize_t x_coord_dim = 8;//totalSize[0];
//                hid_t h5XCorrdSpace = H5Screate_simple(1, &x_coord_dim, 0L);
//                checkH5Err(h5XCorrdSpace);
//                hid_t h5XCoord = H5Dcreate(m_hdfFile,"/x",H5T_IEEE_F64LE, h5XCorrdSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//                checkH5Err(h5XCoord);

//                hsize_t xWriteStrat = 0;
//                hsize_t xWriteCount = 8;

//                hid_t something = H5Screate_simple(2, &xWriteCount, 0L);
//                checkH5Err(something);
//                checkH5Err(H5Sselect_hyperslab(h5XCorrdSpace, H5S_SELECT_SET,
//                        &xWriteStrat, 0L, &xWriteCount, 0L));

//                double x[8] = {0,1,0,1,0,1,0,1};


                //checkH5Err(H5Dwrite(h5XCoord, H5T_NATIVE_DOUBLE, something, h5XCorrdSpace, m_hdfAccessList, x));

                //checkH5Err(H5Sclose(something));
                //checkH5Err(H5Sclose(h5XCoord));
                //checkH5Err(H5Dclose(h5XCorrdSpace));

			// Create partition dataset
			hsize_t partDim = totalSize[0];
			hid_t h5partitionSpace = H5Screate_simple(1, &partDim, 0L);
			checkH5Err(h5partitionSpace);
			hid_t h5partition = H5Dcreate(m_hdfFile, "/partition", H5T_STD_U32LE, h5partitionSpace,
					H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			checkH5Err(h5partition);

			// Create variable datasets
			hsize_t varChunkDims[2] = {
					utils::Env::get<hsize_t>("XDMFWRITER_TIME_CHUNK_SIZE", 1),
					utils::Env::get<hsize_t>("XDMFWRITER_ELEMENT_CHUNK_SIZE", 0)};
			if (varChunkDims[1] == 0)
				// 0 elements -> all elements
				varChunkDims[1] = totalSize[0];
			// TODO add additional check for the chunk size
			hid_t h5pCreate = H5Pcreate(H5P_DATASET_CREATE);
			checkH5Err(h5pCreate);
			checkH5Err(H5Pset_chunk(h5pCreate, 2, varChunkDims));

			for (unsigned int i = 0; i < numVars(); i++) {
				std::string varName("/");
				varName += variables[i];

				hsize_t varDims[2] = {0, totalSize[0]};
				hsize_t varDimsMax[2] = {H5S_UNLIMITED, totalSize[0]};
				hid_t h5VarSpace = H5Screate_simple(2, varDims, varDimsMax);
				checkH5Err(h5VarSpace);

				m_hdfVars[i] = H5Dcreate(m_hdfFile, varName.c_str(), H5T_IEEE_F64LE, h5VarSpace,
						H5P_DEFAULT, h5pCreate, H5P_DEFAULT);
				checkH5Err(m_hdfVars[i]);

				checkH5Err(H5Sclose(h5VarSpace));
			}
			checkH5Err(H5Pclose(h5pCreate));

			// Write connectivity
			hsize_t connectWriteStart[2] = {offset[0], 0};
			hsize_t connectWriteCount[2] = {localSize[0], m_topoTypeSize};
			hid_t h5MemSpace = H5Screate_simple(2, connectWriteCount, 0L);
			checkH5Err(h5MemSpace);
			checkH5Err(H5Sselect_hyperslab(h5connectSpace, H5S_SELECT_SET,
					connectWriteStart, 0L, connectWriteCount, 0L));

			checkH5Err(H5Dwrite(h5connect, H5T_NATIVE_ULONG, h5MemSpace,
					h5connectSpace, m_hdfAccessList, cells));
			checkH5Err(H5Sclose(h5MemSpace));

			// Write geometry
			hsize_t geometryWriteStart[2] = {offset[1], 0};
			hsize_t geometryWriteCount[2] = {localSize[1], 3};
			h5MemSpace = H5Screate_simple(2, geometryWriteCount, 0L);
			checkH5Err(h5MemSpace);
			checkH5Err(H5Sselect_hyperslab(h5geometrySpace, H5S_SELECT_SET,
					geometryWriteStart, 0L, geometryWriteCount, 0L));

			checkH5Err(H5Dwrite(h5geometry, H5T_NATIVE_DOUBLE, h5MemSpace,
					h5geometrySpace, m_hdfAccessList, vertices));
			checkH5Err(H5Sclose(h5MemSpace));

			// Write partition information
			hsize_t partitionWriteStart = offset[0];
			hsize_t partitionWriteCount = localSize[0];
			h5MemSpace = H5Screate_simple(1, &partitionWriteCount, 0L);
			checkH5Err(h5MemSpace);
			checkH5Err(H5Sselect_hyperslab(h5partitionSpace, H5S_SELECT_SET,
					&partitionWriteStart, 0L, &partitionWriteCount, 0L));

			checkH5Err(H5Dwrite(h5partition, H5T_NATIVE_UINT, h5MemSpace,
					h5partitionSpace, m_hdfAccessList, partition));
			checkH5Err(H5Sclose(h5MemSpace));

			// Close all datasets we only write once
			checkH5Err(H5Dclose(h5connect));
			checkH5Err(H5Sclose(h5connectSpace));
			checkH5Err(H5Dclose(h5geometry));
			checkH5Err(H5Sclose(h5geometrySpace));
			checkH5Err(H5Dclose(h5partition));
			checkH5Err(H5Sclose(h5partitionSpace));

			checkH5Err(H5Fflush(m_hdfFile, H5F_SCOPE_GLOBAL));
		} else {
			// Get variable datasets
			for (size_t i = 0; i < numVars(); i++) {
				std::string varName("/");
				varName += variables[i];

				m_hdfVars[i] = H5Dopen(m_hdfFile, varName.c_str(), H5P_DEFAULT);
				checkH5Err(m_hdfVars[i]);
			}
		}

		// Memory space we use for writing one variable at one time step
		hsize_t writeCount[2] = {1, localSize[0]};
		m_hdfVarMemSpace = H5Screate_simple(2, writeCount, 0L);
		checkH5Err(m_hdfVarMemSpace);
	}

	void writeData(unsigned int timestep, unsigned int id, const double* data)
	{
		hsize_t extent[2] = {timestep+1, totalCells()};
		checkH5Err(H5Dset_extent(m_hdfVars[id], extent));

		hid_t h5VarSpace = H5Dget_space(m_hdfVars[id]);
		checkH5Err(h5VarSpace);
		hsize_t writeStart[2] = {timestep, offsetCells()};
		hsize_t writeCount[2] = {1, localCells()};
		checkH5Err(H5Sselect_hyperslab(h5VarSpace, H5S_SELECT_SET, writeStart, 0L, writeCount, 0L));

		checkH5Err(H5Dwrite(m_hdfVars[id], H5T_NATIVE_DOUBLE, m_hdfVarMemSpace,
				h5VarSpace, m_hdfAccessList, data));

		checkH5Err(H5Sclose(h5VarSpace));
	}

	void flush()
	{
		checkH5Err(H5Fflush(m_hdfFile, H5F_SCOPE_GLOBAL));
	}

	/**
	 * Closes the HDF5 file (should be done before MPI_Finalize is called)
	 */
	void close()
	{
		if (m_hdfVars) {
			for (unsigned int i = 0; i < numVars(); i++)
				checkH5Err(H5Dclose(m_hdfVars[i]));
			checkH5Err(H5Sclose(m_hdfVarMemSpace));
			checkH5Err(H5Pclose(m_hdfAccessList));
			checkH5Err(H5Fclose(m_hdfFile));

			delete [] m_hdfVars;
			m_hdfVars = 0L; // Indicates closed file
		}
	}

public:
	static const char* format()
	{
		return "HDF";
	}

	/**
	 * @param prefix
	 * @param var
	 * @return The location name in the XDMF file for a data item
	 */
	static std::string dataItemLocation(const char* prefix, const char* var)
	{
		return std::string(prefix) + fileExention() + ":/" + var;
	}

private:
	/**
	 * @return The file extension for this backend
	 */
	static const char* fileExention()
	{
		return ".h5";
	}

	template<typename T>
	static void checkH5Err(T status)
	{
		if (status < 0)
			logError() << "An HDF5 error occurred in the XDMF writer";
	}
};

}

}

#endif // XDMF_WRITER_BACKENDS_HDF5_H
