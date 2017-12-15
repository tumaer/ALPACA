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

#ifndef XDMF_WRITER_H
#define XDMF_WRITER_H

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>

#ifdef USE_HDF
#endif // USE_HDF

#include "submodules/utils/env.h"

#include "scorep_wrapper.h"
#ifdef PARALLEL
#include "BlockBuffer.h"
#include "ParallelVertexFilter.h"
#else // PARALLEL
#include "BlockBufferSerial.h"
#endif // PARALLEL
#include "backends/Backend.h"

namespace xdmfwriter
{

/**
 * The topology types
 */
enum TopoType {
	TRIANGLE,
    TETRAHEDRON,
    HEXAHEDRON
};

/**
 * Writes data in XDMF format
 */
template<enum TopoType>
class XdmfWriter
{
private:
#ifdef PARALLEL
	MPI_Comm m_comm;
#endif // PARALLEL

	int m_rank;

	std::string m_outputPrefix;

    //std::fstream m_xdmfFile;

	/** The backend for large scale I/O */
	backends::Backend m_backend;

	/** Names of the variables that should be written */
	const std::vector<const char*> m_variableNames;

	/** Total number of cells */
	unsigned long m_totalCells;

	/** Block buffer used to create equal sized blocks */
	BlockBuffer m_blockBuffer;

	/** Buffer required for blocking */
	double *m_blocks;

	/** Offsets in the XDMF file which describe the time dimension size */
	size_t *m_timeDimPos;

	/** Only execute the flush on certain time steps */
	unsigned int m_flushInterval;

	/** Output step counter */
	unsigned int m_timestep;

public:
	/**
	 * @param timestep Set this to > 0 to activate append mode
	 */
	XdmfWriter(int rank, const char* outputPrefix, const std::vector<const char*> &variableNames,
			unsigned int timestep = 0)
		: m_rank(rank), m_outputPrefix(outputPrefix),
		  m_backend(topoTypeSize()),
		  m_variableNames(variableNames),
		  m_totalCells(0),
		  m_blocks(0L), m_timeDimPos(0L), m_flushInterval(0), m_timestep(timestep)
	{
#ifdef PARALLEL
		m_comm = MPI_COMM_WORLD;
#endif // PARALLEL

		if (rank == 0) {
			std::string xdmfName = std::string(outputPrefix) + ".xdmf";

			std::ofstream(xdmfName.c_str(), std::ios::app).close(); // Create the file (if it does not exist)
            //m_xdmfFile.open(xdmfName.c_str());
		}
	}

	virtual ~XdmfWriter()
	{
		if (m_timeDimPos)
			close();
	}

#ifdef PARALLEL
	/**
	 * Sets the communicator that should be used. Default is MPI_COMM_WORLD.
	 */
	void setComm(MPI_Comm comm)
	{
		m_comm = comm;
	}
#endif // PARALLEL

	void init(unsigned int numCells, const unsigned int* cells, unsigned int numVertices, const double *vertices, bool useVertexFilter = true)
	{
        if(useVertexFilter) {;}

		unsigned int offset = 0;
#ifdef PARALLEL
		// Apply vertex filter
		ParallelVertexFilter filter(m_comm);
		if (useVertexFilter) {
			// Filter duplicate vertices
			filter.filter(numVertices, vertices);
			vertices = filter.localVertices();
			numVertices = filter.numLocalVertices();
		} else {
			// No vertex filter -> just get the offset we should at
			offset = numVertices;
			MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, m_comm);
			offset -= numVertices;
		}
#endif // PARALLEL

		// Add vertex offset to all cells and convert to unsigned long
		unsigned long *h5Cells = new unsigned long[numCells * topoTypeSize()];
#ifdef PARALLEL
		if (useVertexFilter) {
#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
#endif
			for (size_t i = 0; i < numCells*topoTypeSize(); i++)
				h5Cells[i] = filter.globalIds()[cells[i]];
		} else
#endif // PARALLEL
		{
#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
#endif
			for (size_t i = 0; i < numCells*topoTypeSize(); i++)
				h5Cells[i] = cells[i] + offset;
		}

		// Initialize the XDMF file
		unsigned long totalSize[2] = {numCells, numVertices};
#ifdef PARALLEL
		MPI_Allreduce(MPI_IN_PLACE, totalSize, 2, MPI_UNSIGNED_LONG, MPI_SUM, m_comm);
#endif // PARALLEL

//		if (m_rank == 0) {
//			std::string backendPrefix(m_outputPrefix);
//			// Remove all directories from prefix
//			size_t pos = backendPrefix.find_last_of('/');
//			if (pos != std::string::npos)
//				backendPrefix = backendPrefix.substr(pos+1);


//			m_xdmfFile << "<?xml version=\"1.0\" ?>" << std::endl
//					<< "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
//					<< "<Xdmf Version=\"2.0\">" << std::endl
//					<< " <Domain>" << std::endl;
//			m_xdmfFile << "  <Topology TopologyType=\"" << topoTypeName() << "\" NumberOfElements=\"" << totalSize[0] << "\">" << std::endl
//					// This should be UInt but for some reason this does not work with binary data
//					<< "   <DataItem NumberType=\"Int\" Precision=\"8\" Format=\""
//						<< backends::Backend::format() << "\" Dimensions=\""
//						<< totalSize[0] << " " << topoTypeSize() << "\">"
//					<< backends::Backend::dataItemLocation(backendPrefix.c_str(), "connect")
//					<< "</DataItem>" << std::endl
//					<< "  </Topology>" << std::endl;
//			m_xdmfFile << "  <Geometry name=\"geo\" GeometryType=\"XYZ\" NumberOfElements=\"" << totalSize[1] << "\">" << std::endl
//					<< "   <DataItem NumberType=\"Float\" Precision=\"8\" Format=\""
//						<< backends::Backend::format() << "\" Dimensions=\"" << totalSize[1] << " 3\">"
//					<< backends::Backend::dataItemLocation(backendPrefix.c_str(), "geometry")
//					<< "</DataItem>" << std::endl
//					<< "  </Geometry>" << std::endl;

//			m_xdmfFile << "  <DataItem NumberType=\"UInt\" Precision=\"4\" Format=\""
//						<< backends::Backend::format() << "\" Dimensions=\"" << totalSize[0] << "\">"
//					<< backends::Backend::dataItemLocation(backendPrefix.c_str(), "partition")
//					<< "</DataItem>" << std::endl;

//			m_timeDimPos = new size_t[m_variableNames.size()];
//			for (size_t i = 0; i < m_variableNames.size(); i++) {
//				m_xdmfFile << "  <DataItem NumberType=\"Float\" Precision=\"8\" Format=\""
//						<< backends::Backend::format() << "\" Dimensions=\"";
//				m_timeDimPos[i] = m_xdmfFile.tellp();
//				m_xdmfFile << std::setw(MAX_TIMESTEP_SPACE) << m_timestep << ' ' << totalSize[0] << "\">"
//						<< backends::Backend::dataItemLocation(backendPrefix.c_str(), m_variableNames[i])
//						<< "</DataItem>" << std::endl;
//			}

//			m_xdmfFile << "  <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;

//			if (m_timestep == 0)
//				closeXdmf();
//			else {
//				// Jump the correct position in the file
//				std::ostringstream tStartStream;
//				timeStepStartXdmf(m_timestep-1, tStartStream);
//				std::string tStart = tStartStream.str();

//				// Find beginning of the (correct) time step
//				std::string line;
//				while (getline(m_xdmfFile, line)) {
//					if (line.find(tStart) != std::string::npos)
//						break;
//				}
//				if (!m_xdmfFile)
//					logError() << "Unable to find time step for appending";

//				// Find end of this time step
//				while (getline(m_xdmfFile, line)) {
//					if (line.find("</Grid>") != std::string::npos)
//						break;
//				}
//			}
//		}

		// Create partition information
		unsigned int *partInfo = new unsigned int[numCells];
#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
		for (unsigned int i = 0; i < numCells; i++)
			partInfo[i] = m_rank;

#ifdef PARALLEL
		// Create block buffer
		unsigned long blockSize = utils::Env::get<unsigned long>("XDMFWRITER_BLOCK_SIZE", 1);
		if (blockSize > 1) {
			m_blockBuffer.init(m_comm, numCells, MPI_DOUBLE, blockSize);

			MPI_Comm newComm;
			MPI_Comm_split(m_comm, m_blockBuffer.count() > 0 ? 1 : MPI_UNDEFINED, 0, &newComm);
			m_comm = newComm; // We no longer need MPI_COMM_WORLD
		}

		// Exchange cells and vertices
		// Be careful with allocation/deallocation!!!
		unsigned long *blockedCells = 0L;
		double *blockedVertices = 0L;
		unsigned int *blockedPartInfo = 0L;
		if (m_blockBuffer.isInitialized()) {
			numCells = m_blockBuffer.count();

			if (m_timestep == 0) {
				// Cells, vertices and partition info only needs to be exchanged
				// when creating a new file
				blockedCells = new unsigned long[numCells * topoTypeSize()];
				m_blockBuffer.exchange(h5Cells, MPI_UNSIGNED_LONG, topoTypeSize(), blockedCells);


				blockedVertices = m_blockBuffer.exchangeAny(vertices, MPI_DOUBLE, 3,
						numVertices, numVertices);

				if (numCells > 0)
					blockedPartInfo = new unsigned int[numCells];
				m_blockBuffer.exchange(partInfo, MPI_UNSIGNED, 1, blockedPartInfo);

				// Overwrite pointers
				delete [] h5Cells;
				h5Cells = blockedCells;

				vertices = blockedVertices;

				delete [] partInfo;
				partInfo = blockedPartInfo;
			}

			// Allocate memory for data exchange
			if (numCells > 0)
				m_blocks = new double[numCells];
		}
#endif // PARALLEL

		// Create and initialize the HDF5 file
		if (m_blockBuffer.count() > 0) {
#ifdef PARALLEL
			m_backend.setComm(m_comm);
#endif // PARALLEL

			// Compute the offsets where we should start in backend file
			unsigned long offsets[2] = {numCells, numVertices};
#ifdef PARALLEL
			MPI_Scan(MPI_IN_PLACE, offsets, 2, MPI_UNSIGNED_LONG, MPI_SUM, m_comm);
#endif // PARALLEL
			offsets[0] -= numCells;
			offsets[1] -= numVertices;

			// Local size
			unsigned int localSize[2] = {numCells, numVertices};

			if (m_timestep == 0) {
				m_backend.open(m_outputPrefix, m_variableNames,
						totalSize, localSize, offsets,
						h5Cells, vertices, partInfo);
#ifdef PARALLEL
				BlockBuffer::free(blockedVertices);
#endif // PARALLEL
				delete [] partInfo;
			} else
				m_backend.open(m_outputPrefix, m_variableNames,
						totalSize, localSize, offsets,
						h5Cells, vertices, partInfo, false);

			// Save total number of cells (required to append a timestep)
			m_totalCells = totalSize[0];

			// Get flush interval
			m_flushInterval = utils::Env::get<unsigned int>("XDMFWRITER_FLUSH_INTERVAL", 1);
		}

		delete [] h5Cells;
	}

	/**
	 * Closes the HDF5 file (should be done before MPI_Finalize is called)
	 */
	void close()
	{
		if (m_blockBuffer.count() > 0) {
			// Close backend
			m_backend.close();

#ifdef PARALLEL
			if (m_blockBuffer.isInitialized())
				MPI_Comm_free(&m_comm);
#endif // PARALLEL
		}

		delete [] m_timeDimPos;
		m_timeDimPos = 0L; // Indicates closed file

		delete [] m_blocks;
	}

	/**
	 * Add a new output time step
	 */
    void addTimeStep(double time) {

    // Avoid Compiler Warnings after butchering this function.
    double a = time;
    a += a;
//		if (m_rank == 0) {
//			m_xdmfFile << "   ";
//			timeStepStartXdmf(m_timestep, m_xdmfFile);
//			m_xdmfFile << std::endl;
//			m_xdmfFile << "    <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>" << std::endl
//					<< "    <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>" << std::endl
//					<< "    <Time Value=\"" << time << "\"/>" << std::endl;
//			m_xdmfFile << "    <Attribute Name=\"partition\" Center=\"Cell\">" << std::endl
//					// Not sure why we need the total cells here but paraview complains otherwise
//					<< "     <DataItem Reference=\"/Xdmf/Domain/DataItem[1]\" Dimensions=\"" << m_totalCells << "\"/>" << std::endl
//					<< "    </Attribute>" << std::endl;
//			for (size_t i = 0; i < m_variableNames.size(); i++) {
//				m_xdmfFile << "    <Attribute Name=\"" << m_variableNames[i] << "\" Center=\"Cell\">" << std::endl
//						<< "     <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << m_totalCells << "\">" << std::endl
//						<< "      <DataItem NumberType=\"UInt\" Precision=\"8\" Format=\"XML\" Dimensions=\"3 2\">"
//						<< m_timestep << " 0 1 1 1 " << m_totalCells << "</DataItem>" << std::endl
//						<< "      <DataItem Reference=\"/Xdmf/Domain/DataItem[" << (i+2) << "]\"/>" << std::endl // PartInfo + 1 based index
//						<< "     </DataItem>" << std::endl
//						<< "    </Attribute>" << std::endl;
//			}
//			m_xdmfFile << "   </Grid>" << std::endl;

//			// Update total steps information
//			size_t pos = m_xdmfFile.tellp();
//			for (size_t i = 0; i < m_variableNames.size(); i++) {
//				m_xdmfFile.seekp(m_timeDimPos[i]);
//				m_xdmfFile << std::setw(MAX_TIMESTEP_SPACE) << (m_timestep+1);
//			}
//			m_xdmfFile.seekp(pos);

//			closeXdmf();
//		}

		m_timestep++;
	}

	/**
	 * Write data for one variable at the current time step
	 *
	 * @param id The number of the variable that should be written
	 */
	void writeData(unsigned int id, const double *data)
	{
		SCOREP_USER_REGION("XDMFWriter_writeData", SCOREP_USER_REGION_TYPE_FUNCTION);

#ifdef PARALLEL
		if (m_blockBuffer.isInitialized()) {
			m_blockBuffer.exchange(data, m_blocks);
			data = m_blocks;
		}
#endif // PARALLEL

		if (m_blockBuffer.count() > 0) {
			m_backend.writeData(m_timestep-1, id, data); // We already incremented m_timestep
		}
	}

	/**
	 * Flushes the data to disk
	 */
	void flush()
	{
		SCOREP_USER_REGION("XDMFWriter_flush", SCOREP_USER_REGION_TYPE_FUNCTION);

		if (m_blockBuffer.count() > 0) {
			if (m_timestep % m_flushInterval == 0)
				m_backend.flush();
		}
	}

	/**
	 * @return The current time step of the output file
	 */
	unsigned int timestep() const
	{
		return m_timestep;
	}

private:
	void closeXdmf()
	{
        //size_t contPos = m_xdmfFile.tellp();
        //m_xdmfFile << "  </Grid>" << std::endl
        //		<< " </Domain>" << std::endl
    //			<< "</Xdmf>" << std::endl;
    //	m_xdmfFile.seekp(contPos);
	}

	/**
	 * @return Name of the topology type in the XDMF file
	 */
	const char* topoTypeName() const;

	/**
	 * @return Number of vertices of the topology type
	 */
	unsigned int topoTypeSize() const;

private:
	/**
	 * Write the beginning of a time step to the stream
	 */
	static void timeStepStartXdmf(unsigned int timestep, std::ostream &s)
	{
		s << "<Grid Name=\"step_" << std::setw(MAX_TIMESTEP_SPACE) << std::setfill('0') << timestep << std::setfill(' ')
				<< "\" GridType=\"Uniform\">";
	}

private:
	static const unsigned int MAX_TIMESTEP_SPACE = 12;
};

template<> inline
const char* XdmfWriter<TRIANGLE>::topoTypeName() const
{
	return "Triangle";
}

template<> inline
const char* XdmfWriter<TETRAHEDRON>::topoTypeName() const
{
	return "Tetrahedron";
}

template<> inline
const char* XdmfWriter<HEXAHEDRON>::topoTypeName() const
{
    return "Hexahedron";
}


template<> inline
unsigned int XdmfWriter<TRIANGLE>::topoTypeSize() const
{
	return 3;
}

template<> inline
unsigned int XdmfWriter<TETRAHEDRON>::topoTypeSize() const
{
	return 4;
}

template<> inline
unsigned int XdmfWriter<HEXAHEDRON>::topoTypeSize() const
{
    return 8;
}

}

#endif // XDMF_WRITER_H
