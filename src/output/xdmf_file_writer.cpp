/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*  \\\\                                                                                  *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA                                                                                 *
* Copyright (c) 2017 Nikolaus A. Adams and contributors (see AUTHORS list)               *
* All rights reserved.                                                                   *
*                                                                                        *
* Chair of Aerodynamics and Fluid Mechanics                                              *
* Technical University of Munich                                                         *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
* This project has received funding from the European Reseach Council (ERC)              *
* under the European Union's Horizon 2020 research and innovation programme              *
* (grant agreement No 667483).                                                           *
*                                                                                        *
* ERC Advanced Grant No 667483, Prof. Dr. Nikolaus A. Adams:                             *
* "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"      *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Redistribution and use in source and binary forms, with or without                     *
* modification, are permitted provided that the following conditions are met:            *
*                                                                                        *
* 1. Redistributions of source code must retain the above copyright notice,              *
*    this list of conditions and the following disclaimer.                               *
*                                                                                        *
* 2. Redistributions in binary form must reproduce the above copyright notice            *
*    this list of conditions and the following disclaimer in the documentation           *
*    and/or other materials provided with the distribution.                              *
*                                                                                        *
* 3. Neither the name of the copyright holder nor the names of its                       *
*    contributors may be used to endorse or promote products derived from this           *
*    software without specific prior written permission.                                 *
*                                                                                        *
* 4. Any redistribution of substantial fractions of the code as a                        *
*    different project should preserve the word ALPACA in the name                       *
*    of the code                                                                         *
*                                                                                        *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"            *
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE              *
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE            *
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE              *
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                    *
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF                   *
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS               *
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN                *
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)                *
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE            *
* POSSIBILITY OF SUCH DAMAGE.                                                            *
*                                                                                        *
* Please note, several third-party tools are used within the ALPACA code under           *
* their own license agreement.                                                           *
*                                                                                        *
* 1. xdmf_writer        : Licensed by Technische Universitaet Muenchen                   *
*                         See 'COPYING_XDMF_WRITER' for more information.                *
*                                                                                        *
* 2. tiny_xml           : This software is provided 'as-is', without any express or      *
*                         implied warranty. In no event will the authors be held         *
*                         liable for any damages arising from the use of this software.  *
*                         See COPYING_TINY_XMLfor more information.                      *
*                                                                                        *
* 3. expression_toolkit : Free use of The C++ Mathematical Expression Toolkit Library is *
*                         permitted under the guidelines and in accordance with the most *
*                         current version of the Common Public License.                  *
*                         http://www.opensource.org/licenses/cpl1.0.php                  *
*                         See COPYING_EXPRESSION_TOOLKITfor more information.            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* AUTHORS                                                                                *
*                                                                                        *
*   Prof. Dr. Nikolaus A. Adams                                                          *
*                                                                                        *
*   Dr. Stefan Adami                                                                     *
*   Vladimir Bogdanov                                                                    *
*   Nico Fleischmann                                                                     *
*   Nils Hoppe                                                                           *
*   Naeimeh Hosseini                                                                     *
*   Jakob Kaiser                                                                         *
*   Aleksandr Lunkov                                                                     *
*   Thomas Paula                                                                         *
*   Josef Winter                                                                         *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
*   nanoshock@aer.mw.tum.de                                                              *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, December 15th 2017                                                             *
*                                                                                        *
*****************************************************************************************/

#include "xdmf_file_writer.h"

#include <fstream>
#include <iostream>
#include <compile_time_constants.h>
#include "mpi.h"

std::vector<double> XdmfFileWriter::output_times_;

/**
 * @brief Default constructor.
 * @param filename The name (including the path) of the xdmf file.
 * @param mpi_rank The MPI rank on which this instance is created.
 * @param print_time_series .
 */
XdmfFileWriter::XdmfFileWriter(const std::string filename, const int mpi_rank, bool print_time_series) :
    series_filename_(filename),
    mpi_rank_(mpi_rank),
    print_time_series_(print_time_series)
{
  if(mpi_rank_ == 0 && print_time_series_) {
    std::ofstream output_stream(series_filename_, std::ios::app);

    //We print the header
    output_stream << "<?xml version=\"1.0\" ?>" << std::endl;
    output_stream << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    output_stream << "<Xdmf Version=\"2.0\">" << std::endl;
    output_stream << " <Domain>" << std::endl;
    output_stream << "  <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;

    output_stream.flush();
    output_stream.close();
  }
}

/**
 * @brief The destructor writes the closing entries into the xdmf file before it destroys itself.
 */
XdmfFileWriter::~XdmfFileWriter() {
  if(mpi_rank_ == 0 && print_time_series_) {
    std::ofstream output_stream(series_filename_, std::ios::app);

    //We print the final clauses
    output_stream << "   <Time TimeType=\"List\">" << std::endl;
    output_stream << "    <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\""<< output_times_.size() <<"\">" << std::endl;
    output_stream << "    ";
    for(const auto& time : output_times_) {
        output_stream << std::to_string(time) << " ";
    }
    output_stream << std::endl;
    output_stream << "    </DataItem>" << std::endl;
    output_stream << "   </Time>" << std::endl;
    output_stream << "  </Grid>" << std::endl;
    output_stream << " </Domain>" << std::endl;
    output_stream << "</Xdmf>" << std::endl;

    output_stream.flush();
    output_stream.close();
  }
}

/**
 * @brief Appends all necessary data for the current time step to the xdmf file for the series representation..
 * @param timestep The current time step.
 * @param hdf5_filename The name of the accompaning HDF5 file (wirtten in another place).
 * @param global_number_of_cells The number of cells
 * @param global_number_of_vertices .
 */
void XdmfFileWriter::AppendTimestep(const std::string hdf5_filename, const unsigned int global_number_of_cells, const unsigned int global_number_of_vertices) const {

    std::string short_file_name = hdf5_filename;
    short_file_name.erase(0,short_file_name.find_last_of("/")+1);

    std::ofstream output_stream(series_filename_, std::ios::app);
    output_stream << "   <Grid Name=\"step_" << output_times_.size() << "\" GridType=\"Uniform\">" << std::endl;
    output_stream << "    <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"" << global_number_of_cells << "\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Int\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << " 8\">" << short_file_name << ":/connect</DataItem>" << std::endl;
    output_stream << "    </Topology>" << std::endl;
    output_stream << "    <Geometry name=\"geometry\" GeometryType=\"XYZ\" NumberOfElements=\"" << global_number_of_vertices << "\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_vertices << " 3\">" << short_file_name << ":/geometry</DataItem>" << std::endl;
    output_stream << "    </Geometry>" << std::endl;
    output_stream << "    <Attribute Name=\"partition\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells <<"\">" << short_file_name << ":/partition</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"density\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/density</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"energy\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/energy</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"x_momentum\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/x_momentum</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    if(CC::DIM() != Dimension::One) {
        output_stream << "    <Attribute Name=\"y_momentum\" Center=\"Cell\">" << std::endl;
        output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/y_momentum</DataItem>" << std::endl;
        output_stream << "    </Attribute>" << std::endl;
    }
    if(CC::DIM() == Dimension::Three) {
        output_stream << "    <Attribute Name=\"z_momentum\" Center=\"Cell\">" << std::endl;
        output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/z_momentum</DataItem>" << std::endl;
        output_stream << "    </Attribute>" << std::endl;
    }
    output_stream << "   </Grid>" << std::endl;

    output_stream.flush();
    output_stream.close();
}

/**
 * @brief Create a a full XDMF file for the current time step, and the current time step only.
 * @param hdf5_filename .
 * @param global_number_of_cells The The number of cells in the output.
 * @param global_number_of_vertices The number of vertices in the output.
 */
void XdmfFileWriter::CreateTimestepFile(const std::string hdf5_filename, const unsigned int global_number_of_cells, const unsigned int global_number_of_vertices) const {

    std::string short_file_name = hdf5_filename;
    short_file_name.erase(0,short_file_name.find_last_of("/")+1);
    std::string xdmf_filename = hdf5_filename.substr(0,hdf5_filename.size()-3) + ".xdmf"; //Cut ".h5" = three caracters.

    std::ofstream output_stream(xdmf_filename, std::ios::trunc);
    output_stream << "<?xml version=\"1.0\" ?>" << std::endl;
    output_stream << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    output_stream << "<Xdmf Version=\"2.0\">" << std::endl;
    output_stream << " <Domain>" << std::endl;
    output_stream << "  <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
    output_stream << "   <Grid Name=\"step_" << output_times_.size() << "\" GridType=\"Uniform\">" << std::endl;
    output_stream << "    <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"" << global_number_of_cells << "\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Int\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << " 8\">" << short_file_name << ":/connect</DataItem>" << std::endl;
    output_stream << "    </Topology>" << std::endl;
    output_stream << "    <Geometry name=\"geometry\" GeometryType=\"XYZ\" NumberOfElements=\"" << global_number_of_vertices << "\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_vertices << " 3\">" << short_file_name << ":/geometry</DataItem>" << std::endl;
    output_stream << "    </Geometry>" << std::endl;
    output_stream << "    <Attribute Name=\"partition\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells <<"\">" << short_file_name << ":/partition</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"density\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/density</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"energy\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/energy</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"x_momentum\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/x_momentum</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    if(CC::DIM() != Dimension::One) {
        output_stream << "    <Attribute Name=\"y_momentum\" Center=\"Cell\">" << std::endl;
        output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/y_momentum</DataItem>" << std::endl;
        output_stream << "    </Attribute>" << std::endl;
    }
    if(CC::DIM() == Dimension::Three) {
        output_stream << "    <Attribute Name=\"z_momentum\" Center=\"Cell\">" << std::endl;
        output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/z_momentum</DataItem>" << std::endl;
        output_stream << "    </Attribute>" << std::endl;
    }
    output_stream << "   </Grid>" << std::endl;
    output_stream << "  </Grid>" << std::endl;
    output_stream << " </Domain>" << std::endl;
    output_stream << "</Xdmf>" << std::endl;

    output_stream.flush();
    output_stream.close();
}

/**
 * @brief Writes a single XDMF file for the current tiem step and appends the series-XDMF with the current step.
 * @param timestep The time step of the current output.
 * @param hdf5_filename The name of the externally (= on another class) created (havy data) h5 file.
 * @param global_number_of_cells The The number of cells in the output.
 * @param global_number_of_vertices The number of vertices in the output.
 */
void XdmfFileWriter::WriteXdmfForTimestep(const double timestep, const std::string hdf5_filename, const unsigned int global_number_of_cells, const unsigned int global_number_of_vertices) const {
  if(mpi_rank_ == 0) {
    output_times_.push_back(timestep);
    AppendTimestep(hdf5_filename,global_number_of_cells,global_number_of_vertices);
    CreateTimestepFile(hdf5_filename,global_number_of_cells,global_number_of_vertices);
  }
}

/**
 * @brief Creates a single XDMF file for the current debug output.
 * @param timestep The time of the current output.
 * @param hdf5_filename .
 * @param global_number_of_cells The number of cells in this output (accross all ranks).
 * @param global_number_of_vertices The number of vertices in this output (accross all ranks).
 */
void XdmfFileWriter::WriteDebugXdmfFile(const double timestep, const std::string hdf5_filename, const unsigned int global_number_of_cells, const unsigned int global_number_of_vertices) const {

  if(mpi_rank_ == 0) {
    std::string short_file_name = hdf5_filename;
    short_file_name.erase(0,short_file_name.find_last_of("/")+1);
    std::string xdmf_filename = hdf5_filename.substr(0,hdf5_filename.size()-3) + ".xdmf"; //Cut ".h5" = three caracters.

    std::ofstream output_stream(xdmf_filename, std::ios::trunc);
    output_stream << "<?xml version=\"1.0\" ?>" << std::endl;
    output_stream << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    output_stream << "<Xdmf Version=\"2.0\">" << std::endl;
    output_stream << " <Domain>" << std::endl;
    output_stream << "  <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
    output_stream << "   <Grid Name=\"step_" << timestep << "\" GridType=\"Uniform\">" << std::endl;
    output_stream << "    <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"" << global_number_of_cells << "\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Int\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << " 8\">" << short_file_name << ":/connect</DataItem>" << std::endl;
    output_stream << "    </Topology>" << std::endl;
    output_stream << "    <Geometry name=\"geometry\" GeometryType=\"XYZ\" NumberOfElements=\"" << global_number_of_vertices << "\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_vertices << " 3\">" << short_file_name << ":/geometry</DataItem>" << std::endl;
    output_stream << "    </Geometry>" << std::endl;
    output_stream << "    <Attribute Name=\"partition\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells <<"\">" << short_file_name << ":/partition</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"density_avg\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/density_old</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"density_rhs\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/density_new</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"energy_avg\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/energy_old</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"energy_rhs\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/energy_new</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"x_momentum_avg\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/x_momentum_old</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    output_stream << "    <Attribute Name=\"x_momentum_rhs\" Center=\"Cell\">" << std::endl;
    output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/x_momentum_new</DataItem>" << std::endl;
    output_stream << "    </Attribute>" << std::endl;
    if(CC::DIM() != Dimension::One) {
        output_stream << "    <Attribute Name=\"y_momentum_avg\" Center=\"Cell\">" << std::endl;
        output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/y_momentum_old</DataItem>" << std::endl;
        output_stream << "    </Attribute>" << std::endl;
        output_stream << "    <Attribute Name=\"y_momentum_rhs\" Center=\"Cell\">" << std::endl;
        output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/y_momentum_new</DataItem>" << std::endl;
        output_stream << "    </Attribute>" << std::endl;
    }
    if(CC::DIM() == Dimension::Three) {
        output_stream << "    <Attribute Name=\"z_momentum_avg\" Center=\"Cell\">" << std::endl;
        output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/z_momentum_old</DataItem>" << std::endl;
        output_stream << "    </Attribute>" << std::endl;
        output_stream << "    <Attribute Name=\"x_momentum_rhs\" Center=\"Cell\">" << std::endl;
        output_stream << "     <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Dimensions=\"" << global_number_of_cells << "\">" << short_file_name << ":/z_momentum_new</DataItem>" << std::endl;
        output_stream << "    </Attribute>" << std::endl;
    }
    output_stream << "   </Grid>" << std::endl;
    output_stream << "  </Grid>" << std::endl;
    output_stream << " </Domain>" << std::endl;
    output_stream << "</Xdmf>" << std::endl;

    output_stream.flush();
    output_stream.close();
  }
}
