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

#include "output_writer.h"

#include <fstream>
#include <iostream>
#include "block.h"
#include "node.h"
#include "compile_time_constants.h"

#include "3rdParty/xdmf_writer/XdmfWriter.h"
#include "3rdParty/xdmf_writer/submodules/utils/logger.h"

/**
 * @brief Creates an object to get the simulation data from the RAM to the hard disk.
 * @param flower The tree to read out the data from.
 * @param topology The topolgy to get information about the MPI status of the simulation.
 * @param setup The setup provides information about user-settings of the simulation.
 * @param io Instance of the I/O manager that handles filesystem access.
 * @param logger The logfile writer for this simulation.
 */
OutputWriter::OutputWriter(const Tree& flower, const TopologyManager& topology, const SimulationSetup& setup, const InputOutputManager& io, LogWriter& logger) :
    tree_(flower),
    topology_(topology),
    setup_(setup),
    input_output_(io),
    xdmf_file_writer_(io.OutputXdmfFileName(),topology_.MyRankId(),XdmfTimeSeries()),
    logger_(logger)
{
}

/**
 * @brief Gives an indication whether or not a timeseries of the outputs in xdmf format should be created.
 * @return True if output format is xdmf, false otherwise.
 */
bool OutputWriter::XdmfTimeSeries() const {
  if(setup_.GetOutputFormat() == OutputType::eXdmf) {
    return true;
  } else {
    return false;
  }
}

/**
 * @brief Gives an increment for the timestamp counter.
 * The increment may be != 1 if the macro timestep covered multiple timestamps at once.
 * @param timestep  The current timestep.
 * @param timestamp_counter The current timestamp counter.
 * @return The increment for the timestamp counter.
 */
unsigned int OutputWriter::TimestampCounterIncrement(const double timestep, const unsigned int timestamp_counter) const {

  const std::vector<double> output_times(setup_.OutputTimes());

  auto first_greater_value = std::upper_bound(output_times.begin(),output_times.end(),timestep);

  return std::distance(output_times.begin(),first_greater_value) - timestamp_counter;

}

/**
 * @brief Gives a decision if at the current time step an output is desired according to the user timestamps
 * @param timestep The current time step.
 * @param timestamp_counter The currently investigated timestamp.
 * @return Decision if output is desired.
 */
bool OutputWriter::OutputTimestep(const double timestep, const unsigned int timestamp_counter) const {

  // In this case no output is necessary
  if(timestep < setup_.OutputTimes()[timestamp_counter]) {
    return false;
  } else {
    return true;
  }
}

/**
 * @brief Triggers the output of the simulation results. Based on user Input the correct type of output is created.
 * @param time_step The time at which the simulation is currently at.
 */
void OutputWriter::WriteOutputFile(const double timestep) const {

  double dimensionalized_time = setup_.DimensionalizeTime(timestep);

  switch (setup_.GetOutputFormat()) {
    case OutputType::eXdmf :
        WriteOutputFileXdmf(dimensionalized_time);
    break;
    case OutputType::eAscii :
        WriteOutputFileAscii(dimensionalized_time);
    break;
    default:
      throw std::invalid_argument("Selected output Type is not possible");
    break;
  }

  if(CC::TEST()) {
    // Rank by rank execution of the non mpi-safe function
    int token = 0;
    if(topology_.MyRankId() == 0) {
        WriteOutputFileAscii_1D(dimensionalized_time);
        if(topology_.NumberOfRanks() > 1) {
            MPI_Send(&token,1,MPI_INT,1,0,MPI_COMM_WORLD);
        }
    } else if(topology_.MyRankId() < topology_.NumberOfRanks()-1) {
        MPI_Recv(&token,1,MPI_INT,topology_.MyRankId()-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        WriteOutputFileAscii_1D(dimensionalized_time);
        MPI_Send(&token,1,MPI_INT,topology_.MyRankId()+1,0,MPI_COMM_WORLD);
    } else {
        MPI_Recv(&token,1,MPI_INT,topology_.MyRankId()-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        WriteOutputFileAscii_1D(dimensionalized_time);
    }
  }

  logger_.LogMessage("Outputfile created at t = " + LogWriter::ConvertDoubleFormat(dimensionalized_time),true,true);

}

/**
 * @brief Creates an output of the current time. For the output two files are created, an xdmf file and either an h5 or a binary file.
 *        This setting is made in main.cpp. The created file is named according to the time step. Only writes conservative values.
 *        $Only considers Average Buffers, therefore it may only be called afer a time step has been completed on level zero$.
 * @param timestep The current time of the simulation.
 */
void OutputWriter::WriteOutputFileXdmf(const double timestep) const {

  unsigned int number_of_cells;
  unsigned int number_of_vertices;

  std::array<double, 3> coordinates;

  std::vector<const char*> variable_names;
  variable_names.push_back("density");
  variable_names.push_back("energy");
  variable_names.push_back("x_momentum");
  if(CC::DIM() != Dimension::One) {
    variable_names.push_back("y_momentum");
  }
  if(CC::DIM() == Dimension::Three) {
    variable_names.push_back("z_momentum");
  }

  double block_size;

  std::string output_filename = input_output_.OutputDomainBaseName() +"/data_" + std::to_string(timestep*setup_.GetTimeNamingFactor());

  number_of_cells = tree_.Leaves().size() * CC::ICX() * CC::ICY() * CC::ICZ();
  number_of_vertices = tree_.Leaves().size() * (CC::ICX()+1) * (CC::ICY()+1) * (CC::ICZ()+1);

  //NH Vectors needed here, otherwise stack overflow possible. Vector automatically goes onto the heap.
  std::vector<std::vector<double>> data(variable_names.size(), std::vector<double> (number_of_cells, 0.0));

  std::vector<unsigned int> cells(number_of_cells * 8); // Eight Points define a cube in 3D
  std::vector<double> vertices(number_of_vertices * 3); // 3 as for XYZ

  unsigned int leaves_counter = 0;
  unsigned int vertices_coordinates_counter = 0;
  unsigned int vertices_counter = 0;

  // this loop is needed to know apriori the number of leaves on the current rank
  for(const auto& node : tree_.Leaves()) {
    const Block& block = node->GetBlock();
    block_size = setup_.DimensionalizeLength(node->GetBlockSize());
    coordinates = node->GetBlockDomainCoordinates();
    coordinates[0] = setup_.DimensionalizeLength(coordinates[0]);
    coordinates[1] = setup_.DimensionalizeLength(coordinates[1]);
    coordinates[2] = setup_.DimensionalizeLength(coordinates[2]);

    for(unsigned int k = 0; k <= CC::ICZ(); ++k) {
      for(unsigned int j = 0; j <= CC::ICY(); ++j ) {
        for(unsigned int i = 0; i <= CC::ICX(); ++i) {
          vertices[ (vertices_coordinates_counter + 0)] = coordinates[0] + (double(i)/double(CC::ICX())) * block_size;
          vertices[ (vertices_coordinates_counter + 1)] = coordinates[1] + (double(j)/double(CC::ICY())) * block_size;
          vertices[ (vertices_coordinates_counter + 2)] = coordinates[2] + (double(k)/double(CC::ICZ())) * block_size;
          // update the counter skipping the 3 elements
          vertices_coordinates_counter+=3;
        }
      }
    }
    // We do not want to include the last vertex in each dimension, because it does not build a cell
    for(unsigned int k = 0; k < CC::ICZ(); ++k) {
      for(unsigned int j = 0; j < CC::ICY(); ++j ) {
        for(unsigned int i = 0; i < CC::ICX(); ++i) {
          cells[(vertices_counter + 0)] = i 	+ j 	* (CC::ICX()+1) + k	 	* (CC::ICX()+1) * (CC::ICY()+1) + leaves_counter * (CC::ICX()+1) * (CC::ICY()+1) * (CC::ICZ()+1);
          cells[(vertices_counter + 1)] = (i+1) + j 	* (CC::ICX()+1) + k 	* (CC::ICX()+1) * (CC::ICY()+1) + leaves_counter * (CC::ICX()+1) * (CC::ICY()+1) * (CC::ICZ()+1);
          cells[(vertices_counter + 2)] = (i+1) + (j+1) * (CC::ICX()+1) + k 	* (CC::ICX()+1) * (CC::ICY()+1) + leaves_counter * (CC::ICX()+1) * (CC::ICY()+1) * (CC::ICZ()+1);
          cells[(vertices_counter + 3)] = i 	+ (j+1) * (CC::ICX()+1) + k 	* (CC::ICX()+1) * (CC::ICY()+1) + leaves_counter * (CC::ICX()+1) * (CC::ICY()+1) * (CC::ICZ()+1);
          cells[(vertices_counter + 4)] = i 	+ j 	* (CC::ICX()+1) + (k+1) * (CC::ICX()+1) * (CC::ICY()+1) + leaves_counter * (CC::ICX()+1) * (CC::ICY()+1) * (CC::ICZ()+1);
          cells[(vertices_counter + 5)] = (i+1) + j 	* (CC::ICX()+1) + (k+1) * (CC::ICX()+1) * (CC::ICY()+1) + leaves_counter * (CC::ICX()+1) * (CC::ICY()+1) * (CC::ICZ()+1);
          cells[(vertices_counter + 6)] = (i+1) + (j+1) * (CC::ICX()+1) + (k+1) * (CC::ICX()+1) * (CC::ICY()+1) + leaves_counter * (CC::ICX()+1) * (CC::ICY()+1) * (CC::ICZ()+1);
          cells[(vertices_counter + 7)] = i 	+ (j+1) * (CC::ICX()+1) + (k+1) * (CC::ICX()+1) * (CC::ICY()+1) + leaves_counter * (CC::ICX()+1) * (CC::ICY()+1) * (CC::ICZ()+1);

		  vertices_counter += 8;
        }
      }
    }

    const double (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(0);
    for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                data[0][ (i - CC::FICX())
                       + (j - CC::FICY()) * CC::ICX()
                       + (k - CC::FICZ()) * CC::ICX() * CC::ICY()
                       + leaves_counter * CC::ICX() * CC::ICY() * CC::ICZ()]
                       = setup_.DimensionalizeRho(density[i][j][k]);
            }
        }
    }

    const double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(1);
    for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                data[1][ (i - CC::FICX())
                       + (j - CC::FICY()) * CC::ICX()
                       + (k - CC::FICZ()) * CC::ICX() * CC::ICY()
                       + leaves_counter * CC::ICX() * CC::ICY() * CC::ICZ()]
                       = setup_.DimensionalizeEnergy(energy[i][j][k]);
            }
        }
    }

    for(unsigned int e = 2; e < CC::NoEq(); ++e){
      const double (&cell_data)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(e);
      for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
          for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
            data[e][ (i - CC::FICX())
                   + (j - CC::FICY()) * CC::ICX()
                   + (k - CC::FICZ()) * CC::ICX() * CC::ICY()
                   + leaves_counter * CC::ICX() * CC::ICY() * CC::ICZ()]
                   = setup_.DimensionalizeMomentum(cell_data[i][j][k]);
          }
        }
      }
    }
    leaves_counter++;
  }

  xdmfwriter::XdmfWriter<xdmfwriter::HEXAHEDRON> writer(int(topology_.MyRankId()), output_filename.c_str(), variable_names, 0);

  writer.init(number_of_cells, cells.data(), number_of_vertices, vertices.data(), false);

  writer.addTimeStep(timestep);

  for(unsigned int i = 0; i < variable_names.size(); ++i){
    writer.writeData(i, data[i].data());
  }

  writer.flush();
  writer.close();

  unsigned int global_number_of_cells = 0;
  unsigned int global_number_of_vertices = 0;

  MPI_Reduce(&number_of_cells, &global_number_of_cells,1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&number_of_vertices, &global_number_of_vertices, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  xdmf_file_writer_.WriteXdmfForTimestep(timestep,output_filename+".h5",global_number_of_cells,global_number_of_vertices);

}

/**
 * @brief Creates an output of the current time. The output format is in ASCII for Tecplot. The created file is named according to the
 *        time step. Only writes conservative values. $Only considers Avergae Buffers, therefore it may only be called afer a time step has been completed on level zero$.
 * @param time_step The current time of the simulation.
 */
void OutputWriter::WriteOutputFileAscii(const double timestep) const {

  double block_size;
  std::array<double, 3> coordinates;


  std::string filename = input_output_.OutputDomainBaseName() + "/domain_" + std::to_string(topology_.MyRankId()) + "_"
                       + std::to_string(timestep*setup_.GetTimeNamingFactor()) + ".dat";

  std::ofstream output_stream(filename, std::ios::app);

  for(const auto& node : tree_.Leaves()) {

    coordinates = node->GetBlockDomainCoordinates();
    coordinates[0] = setup_.DimensionalizeLength(coordinates[0]);
    coordinates[1] = setup_.DimensionalizeLength(coordinates[1]);
    coordinates[2] = setup_.DimensionalizeLength(coordinates[2]);
    block_size = setup_.DimensionalizeLength(node->GetBlockSize());


    const Block& block = node->GetBlock();
    unsigned int output_number_x = 9;
    unsigned int output_number_y = 9;
    unsigned int output_number_z = 9;
    std::string momentum_in_x = "rhoU ";
    std::string momentum_in_y = "rhoV ";
    std::string momentum_in_z = "rhoW ";
    if (CC::DIM() == Dimension::One){
        output_number_y = 2;
        output_number_z = 2;
        momentum_in_y = "";
        momentum_in_z = "";

    }
    if (CC::DIM() == Dimension::Two){
        output_number_z = 2;
        momentum_in_z = "";
    }
    output_stream << "title='View' \n";
    output_stream << "variables=x, y, z, rho E "<<momentum_in_x<<momentum_in_y<<momentum_in_z<<"rank \n";
    output_stream << "zone t='" << node->GetId() << "' i= "<<output_number_x<<"  j= "<<output_number_y<<" k= "<<output_number_z;
    output_stream <<" DATAPACKING=BLOCK, VARLOCATION=([4-9]=CELLCENTERED), SOLUTIONTIME = "<< timestep << "\n";

    //NH: We need one more point than cells in a direction
    for(unsigned int k = 0; k < CC::ICZ()+1; ++k) {
        for(unsigned int j = 0; j < CC::ICY()+1; ++j ) {
            for(unsigned int i = 0; i < CC::ICX()+1; ++i) {
                output_stream << coordinates[0] + (double(i)/double(CC::ICX())) * block_size << " ";
            }
            output_stream << "\n";
        }
    }

    output_stream << "\n";

    for(unsigned int k = 0; k < CC::ICZ()+1 ; ++k) {
        for(unsigned int j = 0; j < CC::ICY()+1; ++j ) {
            for(unsigned int i = 0; i < CC::ICX()+1; ++i) {
                output_stream << coordinates[1] + (double(j)/double(CC::ICY())) * block_size << " ";
            }
            output_stream << "\n";
        }
    }

    output_stream << "\n";

    for(unsigned int k = 0; k < CC::ICZ()+1; ++k) {
        for(unsigned int j = 0; j < CC::ICY()+1; ++j ) {
            for(unsigned int i = 0; i < CC::ICX()+1; ++i) {
                output_stream << coordinates[2] + (double(k)/double(CC::ICZ())) * block_size << " ";
            }
            output_stream << "\n";
        }
    }

    output_stream << "\n";

    const double (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(0);
    for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                output_stream << setup_.DimensionalizeRho(density[i][j][k]) << " ";
            }
            output_stream << "\n";
        }
    }

    output_stream << "\n";

    const double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(1);
    for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                output_stream << setup_.DimensionalizeEnergy(energy[i][j][k]) << " ";
            }
            output_stream << "\n";
        }
    }

    output_stream << "\n";

    for(int e = 2; e < CC::NoEq(); ++e) {
        const double (&momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(e);
        for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
                for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                    output_stream << setup_.DimensionalizeMomentum(momentum[i][j][k]) << " ";
                }
                output_stream << "\n";
            }
        }
    output_stream << "\n";
    }

    for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                output_stream << topology_.GetRankOfNode(node->GetId()) << " ";
            }
            output_stream << "\n";
        }
    }

  }

  output_stream.flush();
  output_stream.close();
}


/**
 * @brief Creates an output of the current time for quasi-1D cases in x-direction. The output format is in ASCII for Tecplot.
 *        The created file is named according to the time step. Only writes conservative values. $NOT MPI SAFE, MAY ONLY BE
 *        CALLED BY ONE RANK AT A TIME$ $Only considers Average Buffers, therefore it may only be called afer a time step has
 *        been completed on level zero$.
 * @param time_step The current time of the simulation.
 */
void OutputWriter::WriteOutputFileAscii_1D(const double timestep) const {

  double cell_size;
  unsigned int level;
  std::array<double, 3> coordinates;

  std::string filename = input_output_.OutputFolderName() + "/1DOutput/1d_output_ALTS_" + std::to_string(timestep*setup_.GetTimeNamingFactor()) + ".dat";

  std::ofstream output_stream(filename, std::ios::app);

  for(const auto& node : tree_.Leaves()) {
    coordinates = node->GetBlockDomainCoordinates();
    cell_size = node->GetCellSize();
    level     = node->GetLevel();
    if(coordinates[1] == 0 && coordinates[2]==0) {
        const Block& block = node->GetBlock();
        unsigned int j = CC::FICY();
        unsigned int k = CC::FICZ();
        for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
            output_stream << coordinates[0] + (double(i-4)+0.5) * cell_size << ", ";
            output_stream << block.GetAverageBuffer(0)[i][j][k] << ", ";
            output_stream << block.GetAverageBuffer(1)[i][j][k] << ", ";
            output_stream << block.GetAverageBuffer(2)[i][j][k] << ", ";
            output_stream << block.GetAverageBuffer(3)[i][j][k] << ", ";
            output_stream << block.GetAverageBuffer(4)[i][j][k] << ", ";
            output_stream << cell_size                          << ", ";
            output_stream << level                              << "\n";
        }

    }

  }

  output_stream.flush();
  output_stream.close();
}

/**
 * @brief Creates a restart file at the current time of the simulation, which allows to re-initialize the current state and use as
 *        starting point for a new simulation. $Only considers Average Buffers, therefore it may only be called afer a time step has been completed on level zero$.
 *        %CURRENTLY RESTART OPTION IS NOT AVAILABLE%
 * @param time_step The current time of the simulation.
 */
void OutputWriter::WriteRestartFile(const double timestep) const {

    /* NH 2017-02-02: Currently Function is disaled /still a stub from old tecplot output. Will change in the future.
   * Stub left over as it needs only minor modifications to become a proper restart writer. Less work required n keeping than in re-writing
   */
    throw std::logic_error("Restart Function not yet wroking");

    //FILENAME INCOORECT!!!
    std::string filename = input_output_.OutputFolderName() + "/leveled_ghosts/allGhosts_" + std::to_string(topology_.MyRankId()) + "_"
                         + LogWriter::ConvertDoubleFormat(timestep*setup_.GetTimeNamingFactor()) + ".rst";

    std::ofstream output_stream(filename, std::ios::app);

    for(const auto& level : tree_.FullNodeList() ) {
        for(const auto& node : level) {

            const Block& block = node->GetBlock();

            //NH This tecplot stuff hast ot got, Node id and some other Info however hast to stay.
            output_stream << "title='View' \n";
            output_stream << "variables=x, y, z, rho, E, rhoU, rhoV, rhoW, rank \n";
            output_stream << "zone t='" << node->GetId() << "' i= 17  j= 17 k= 17 DATAPACKING=BLOCK, VARLOCATION=([4-9]=CELLCENTERED), SOLUTIONTIME = "<< timestep << "\n";

            for(int e = 0; e < CC::NoEq(); ++e) {
                const double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(e);
                for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                    for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                        for(unsigned int i = 0; i < CC::TCX(); ++i) {
                            output_stream << cells[i][j][k] << " ";
                        }
                        output_stream << "\n";
                    }
                }
            }

            output_stream << "\n";

            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int i = 0; i < CC::TCX(); ++i) {
                        output_stream << topology_.GetRankOfNode(node->GetId()) << " ";
                    }
                    output_stream << "\n";
                }
            }

        }
    }
}

/**
 * @brief Triggers the creation of an outputfile in the format specified in the inputfile.
 * @param time The 'time' of this debug output. Commonly only an integer number.
 */
void OutputWriter::WriteDebugFile(const double time) const {

  switch (setup_.GetOutputFormat()) {
    case OutputType::eXdmf :
        WriteDebugXdmf(time);
    break;
    case OutputType::eAscii :
        WriteDebugAscii(time);
    break;
    default:
      throw std::invalid_argument("Selected output Type is not possible");
    break;
  }

   logger_.LogMessage("Debug file debug_" + LogWriter::ConvertDoubleFormat(time) + " written");
}

/**
 * @brief Creates an output at the current state of the simulation. The output format is in ASCII for Tecplot. Writes out all available data,
 *        therefore it can be used at any time in the simulation. However, files are very huge and should not be used extensively and for debugging purposes only.
 * @param time 'Time' of the debug output. Is not dimensionalized or modified in any way!
 */
void OutputWriter::WriteDebugAscii(const double time) const {

  double block_size;
  std::array<double, 3> coordinates;
  double z_offset = 0;

  std::string filename = input_output_.OutputDebugBaseName() + "/debug_" + std::to_string(topology_.MyRankId()) + "_"
                       + LogWriter::ConvertDoubleFormat(time) + ".dat";

  std::ofstream output_stream(filename, std::ios::app);

  double dimension_offset_factor = 1.0;
  if(CC::DIM() != Dimension::Three) {dimension_offset_factor = 0.5;}
  unsigned int output_number_x = CC::TCX()+1;
  unsigned int output_number_y = CC::TCY()+1;
  unsigned int output_number_z = CC::TCZ()+1;
  std::string momentum_in_x = "rhoU";
  std::string momentum_in_y = "";
  std::string momentum_in_z = "";
  std::string new_momentum_in_x = "rhoU_new";
  std::string new_momentum_in_y = "";
  std::string new_momentum_in_z = "";

  if (CC::DIM() != Dimension::One){
    momentum_in_y = "rhoV";
    new_momentum_in_y = "rhoV_new";
  }

  if (CC::DIM() == Dimension::Three){
    momentum_in_z = "rhoW";
    new_momentum_in_z = "rhoW_new";
  }

  for(const auto& level_iterator : tree_.FullNodeList() ) {
     for(const auto& it : level_iterator) {
        coordinates = it->GetBlockDomainCoordinates();
        coordinates[0] = setup_.DimensionalizeLength(coordinates[0]);
        coordinates[1] = setup_.DimensionalizeLength(coordinates[1]);
        coordinates[2] = setup_.DimensionalizeLength(coordinates[2]);
        block_size = setup_.DimensionalizeLength(it->GetBlockSize());

        const Block& block = it->GetBlock();

        output_stream << "variables=x, y, z, rho, E, "<< momentum_in_x << " " << momentum_in_y << " " << momentum_in_z << " rho_new, E_new, " << new_momentum_in_x;
        output_stream << " " << new_momentum_in_y << " " << new_momentum_in_z << " rank \n";
        output_stream << "zone t='" << it->GetId() << "' i= "<< output_number_x << " j= " << " " << output_number_y << " k= " << output_number_z;
        output_stream << " DATAPACKING=BLOCK, VARLOCATION=([4-14]=CELLCENTERED), SOLUTIONTIME = ";
        output_stream << std::scientific << std::setprecision(54) << time << "\n";

        // NH We need one more point than cells
        for(unsigned int k = 0; k < CC::TCZ()+1; ++k) {
            for(unsigned int j = 0; j < CC::TCY()+1; ++j ) {
                for(unsigned int i = 0; i < CC::TCX()+1; ++i) {
                    output_stream << std::scientific << std::setprecision(54) << coordinates[0] + (double(i)/double(CC::TCX()+1)) * block_size << " ";
                }
                output_stream << "\n";
            }
        }

        output_stream << "\n";

        for(unsigned int k = 0; k < CC::TCZ()+1; ++k) {
            for(unsigned int j = 0; j < CC::TCY()+1; ++j ) {
                for(unsigned int i = 0; i < CC::TCX()+1; ++i) {
                    output_stream << std::scientific << std::setprecision(54) << coordinates[1] + (double(j)/double(CC::TCY()+1)) * block_size << " ";
                }
                output_stream << "\n";
            }
        }

        output_stream << "\n";

        for(unsigned int k = 0; k < CC::TCZ()+1; ++k) {
            for(unsigned int j = 0; j < CC::TCY()+1; ++j ) {
                for(unsigned int i = 0; i < CC::TCX()+1; ++i) {
                    output_stream << std::scientific << std::setprecision(54) << z_offset + coordinates[2] + (double(k)/double(CC::TCZ()+1)) * block_size << " ";
                }
                output_stream << "\n";
            }
        }

        output_stream << "\n";

        const double (&density_average_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(0);
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int i = 0; i < CC::TCX(); ++i) {
                        output_stream << std::scientific << std::setprecision(54) << setup_.DimensionalizeRho(density_average_buffer[i][j][k]) << " ";
                    }
                output_stream << "\n";
            }
        }

        output_stream << "\n";

        const double (&energy_average_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(1);
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int i = 0; i < CC::TCX(); ++i) {
                        output_stream << std::scientific << std::setprecision(54) << setup_.DimensionalizeEnergy(energy_average_buffer[i][j][k]) << " ";
                    }
                output_stream << "\n";
            }
        }

        output_stream << "\n";

        for(int e = 2; e < CC::NoEq(); ++e) {
            const double (&momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(e);
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int i = 0; i < CC::TCX(); ++i) {
                        output_stream << std::scientific << std::setprecision(54) << setup_.DimensionalizeMomentum(momentum[i][j][k]) << " ";
                    }
                    output_stream << "\n";
                }
            }
            output_stream << "\n";
        }

        output_stream << "\n";

        const double (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(0);
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int i = 0; i < CC::TCX(); ++i) {
                        output_stream << std::scientific << std::setprecision(54) << setup_.DimensionalizeRho(density[i][j][k]) << " ";
                    }
                output_stream << "\n";
            }
        }

        output_stream << "\n";

        const double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(1);
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int i = 0; i < CC::TCX(); ++i) {
                        output_stream << std::scientific << std::setprecision(54) << setup_.DimensionalizeEnergy(energy[i][j][k]) << " ";
                    }
                output_stream << "\n";
            }
        }

        output_stream << "\n";

        for(int e = 2; e < CC::NoEq(); ++e) {
            const double (&momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(e);
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int i = 0; i < CC::TCX(); ++i) {
                        output_stream << std::scientific << std::setprecision(54) << setup_.DimensionalizeMomentum(momentum[i][j][k]) << " ";
                    }
                    output_stream << "\n";
                }
            }
            output_stream << "\n";
        }

        output_stream << "\n";

        for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                for(unsigned int i = 0; i < CC::TCX(); ++i) {
                    output_stream << std::scientific << std::setprecision(54) << topology_.GetRankOfNode(it->GetId()) << " ";
                }
                output_stream << "\n";
            }
        }

      }
      z_offset += (setup_.GetLevelZeroBlocksZ() + 1) * (setup_.LevelZeroBlockSize()) * dimension_offset_factor;
   }

   output_stream.flush();
   output_stream.close();
}

/**
 * @brief Writes the used time steps into a separate .txt file.
 * @param timesteps_on_finest_level Time steps used on the finest level to be saved in the file.
 */
// NH 2017-02-20 May go out in later Version!
void OutputWriter::WriteTimestepFile(const std::vector<double> timesteps_on_finest_level) const {

  if(topology_.MyRankId() == 0) {
    std::string filename = input_output_.OutputFolderName() + "/time.txt";
    std::ofstream output_stream(filename,std::ios::app);
     for(const auto& timestep : timesteps_on_finest_level) {
        output_stream << timestep << "\n";
     }

    output_stream.flush();
    output_stream.close();

  }
}

/**
 * @brief Creates a debug file in xdmf/h5 format.
 * @param time 'Time' of the debug output. Is not dimensionalized or modified in any way!
 */
void OutputWriter::WriteDebugXdmf(const double time) const {

  unsigned int number_of_cells, number_of_vertices;

  std::array<double, 3> coordinates;

  std::vector<const char*> variables;
  variables.push_back("density_old");
  variables.push_back("density_new");
  variables.push_back("energy_old");
  variables.push_back("energy_new");
  variables.push_back("x_momentum_old");
  variables.push_back("x_momentum_new");
  if(CC::DIM() != Dimension::One) {
    variables.push_back("y_momentum_old");
    variables.push_back("y_momentum_new");
  }
  if(CC::DIM() == Dimension::Three) {
    variables.push_back("z_momentum_old");
    variables.push_back("z_momentum_new");
  }

  double block_size;
  double z_offset = 0;
  double dimension_offset_factor = 1.0;
  if(CC::DIM() != Dimension::Three) {dimension_offset_factor = 0.5;}
  unsigned int number_of_blocks = 0;

  for(const auto& level: tree_.FullNodeList()) {
    number_of_blocks += level.size();
  }

  std::string filename = input_output_.OutputDebugBaseName() + "/debug_" + std::to_string(time);

  number_of_cells = number_of_blocks * CC::TCX() * CC::TCY() * CC::TCZ();
  number_of_vertices = number_of_blocks * (CC::TCX()+1) * (CC::TCY()+1) * (CC::TCZ()+1);

  std::vector< std::vector<double> > data(variables.size(), std::vector<double> (number_of_cells, 0.0));

  std::vector<unsigned int> point_connectivity(number_of_cells * 8);
  std::vector<double> vertices(number_of_vertices * 3);

  unsigned int leaves_counter = 0;
  unsigned int vertices_coord_counter = 0;
  unsigned int vertices_counter = 0;

  static constexpr unsigned int offset = (CC::TCX()+1) * (CC::TCY()+1) * (CC::TCZ()+1);

  for(const auto& level_iterator : tree_.FullNodeList() ) {
    for(const auto& it : level_iterator) {
      const Block block = it->GetBlock();
      block_size = it->GetBlockSize();
      coordinates = it->GetBlockDomainCoordinates();

      for(unsigned int k = 0; k <= CC::TCZ(); ++k) {
        for(unsigned int j = 0; j <= CC::TCY(); ++j ) {
          for(unsigned int i = 0; i <= CC::TCX(); ++i) {
            vertices[ (vertices_coord_counter + 0)] =            coordinates[0] + (double(i)/double(CC::TCX()+1)) * block_size;
            vertices[ (vertices_coord_counter + 1)] =            coordinates[1] + (double(j)/double(CC::TCY()+1)) * block_size;
            vertices[ (vertices_coord_counter + 2)] = z_offset + coordinates[2] + (double(k)/double(CC::TCZ()+1)) * block_size;
            // update the counter skipping the 3 elements
            vertices_coord_counter+=3;
          }
        }
      }

      // We do not want to include the last vertex in each dimension, because it does not build a cell
      for(unsigned int k = 0; k < CC::TCZ(); ++k) {
        for(unsigned int j = 0; j < CC::TCY(); ++j ) {
          for(unsigned int i = 0; i < CC::TCX(); ++i) {
            point_connectivity[(vertices_counter + 0)] = i     + j     * (CC::TCX()+1) + k	   * (CC::TCX()+1) * (CC::TCY()+1) + leaves_counter * offset;
            point_connectivity[(vertices_counter + 1)] = (i+1) + j     * (CC::TCX()+1) + k 	   * (CC::TCX()+1) * (CC::TCY()+1) + leaves_counter * offset;
            point_connectivity[(vertices_counter + 2)] = (i+1) + (j+1) * (CC::TCX()+1) + k 	   * (CC::TCX()+1) * (CC::TCY()+1) + leaves_counter * offset;
            point_connectivity[(vertices_counter + 3)] = i     + (j+1) * (CC::TCX()+1) + k 	   * (CC::TCX()+1) * (CC::TCY()+1) + leaves_counter * offset;
            point_connectivity[(vertices_counter + 4)] = i     + j     * (CC::TCX()+1) + (k+1) * (CC::TCX()+1) * (CC::TCY()+1) + leaves_counter * offset;
            point_connectivity[(vertices_counter + 5)] = (i+1) + j     * (CC::TCX()+1) + (k+1) * (CC::TCX()+1) * (CC::TCY()+1) + leaves_counter * offset;
            point_connectivity[(vertices_counter + 6)] = (i+1) + (j+1) * (CC::TCX()+1) + (k+1) * (CC::TCX()+1) * (CC::TCY()+1) + leaves_counter * offset;
            point_connectivity[(vertices_counter + 7)] = i     + (j+1) * (CC::TCX()+1) + (k+1) * (CC::TCX()+1) * (CC::TCY()+1) + leaves_counter * offset;

            vertices_counter += 8;
          }
        }
      }

      for(unsigned int e = 0; e < CC::NoEq(); ++e){
        // old equations
        const double (&cell_data_old)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(e);
        for(unsigned int k = 0; k < CC::TCZ(); ++k) {
          for(unsigned int j = 0; j < CC::TCY(); ++j ) {
            for(unsigned int i = 0; i < CC::TCX(); ++i) {
                data[e*2][i + j * CC::TCX() + k * CC::TCX() * CC::TCY() + leaves_counter * CC::TCX() * CC::TCY() * CC::TCZ()] = cell_data_old[i][j][k];
            }
          }
        }

        // new equations
        const double (&cell_data_new)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(e);
        for(unsigned int k = 0; k < CC::TCZ(); ++k) {
          for(unsigned int j = 0; j < CC::TCY(); ++j ) {
            for(unsigned int i = 0; i < CC::TCX(); ++i) {
              data[e*2 + 1][i + j * CC::TCX() + k * CC::TCX() * CC::TCY() + leaves_counter * CC::TCX() * CC::TCY() * CC::TCZ()] = cell_data_new[i][j][k];
            }
          }
        }
      }

      leaves_counter++;
    }
    z_offset += (setup_.GetLevelZeroBlocksZ() + 1) * (setup_.LevelZeroBlockSize()) * dimension_offset_factor;
  }

  xdmfwriter::XdmfWriter<xdmfwriter::HEXAHEDRON> writer(int(topology_.MyRankId()), filename.c_str(), variables, 0);

  writer.init(number_of_cells, point_connectivity.data(), number_of_vertices, vertices.data(), false);
  writer.addTimeStep(time);

  for(unsigned int i = 0; i < variables.size(); ++i){
    writer.writeData(i, data[i].data());
  }

  writer.flush();
  writer.close();

  unsigned int global_number_of_cells = 0;
  unsigned int global_number_of_vertices = 0;

  MPI_Reduce(&number_of_cells, &global_number_of_cells,1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&number_of_vertices, &global_number_of_vertices, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  xdmf_file_writer_.WriteDebugXdmfFile(time*setup_.GetTimeNamingFactor(),filename+".h5",global_number_of_cells,global_number_of_vertices);
}
