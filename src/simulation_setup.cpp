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

#include "simulation_setup.h"

#include <numeric>
#include "mpi.h"
#include "materials/material_manager.h"

/**
 * @brief Default constructor to create a simulation setup from an input file.
 * @param io The I/O manager used to read the inputfile.
 * @param logger Reference to a logger instance to save important messages.
 */
SimulationSetup::SimulationSetup(const InputOutputManager& io, LogWriter& logger) :
    input_output_(io),
    initial_condition_(io.GetInputFileParser()),
    unit_handler_(io.GetInputFileParser(),logger),
    output_format_(io.GetInputFileParser().ReadOutputFormat()),
    output_times_(DetermineOutputTimes(io.GetInputFileParser())),  
    start_time_(unit_handler_.NonDimensionalizeTime(io.GetInputFileParser().ReadStartTime())),
    end_time_(unit_handler_.NonDimensionalizeTime(io.GetInputFileParser().ReadEndTime())),
    time_naming_factor_(io.GetInputFileParser().ReadTimeNamingFactor()),
    maximum_level_(io.GetInputFileParser().ReadMaximumLevel()),
    x_number_of_level_zero_blocks_(io.GetInputFileParser().ReadNumberOfBlocksX()),
    y_number_of_level_zero_blocks_(CC::DIM() != Dimension::One   ? io.GetInputFileParser().ReadNumberOfBlocksY() : 1),
    z_number_of_level_zero_blocks_(CC::DIM() == Dimension::Three ? io.GetInputFileParser().ReadNumberOfBlocksZ() : 1),
    block_size_on_level_zero_(unit_handler_.NonDimensionalizeLength(io.GetInputFileParser().ReadBlockSize())),
    boundary_conditions_(io.GetInputFileParser().ReadBoundaryConditions()),
    epsilon_reference_(ComputeReferenceEpsilon(io.GetInputFileParser().ReadReferenceMultiResolutionEpsilon(),io.GetInputFileParser().ReadLevelOfReferenceMultiResolutionEpsilon())),
    gravity_(unit_handler_.NonDimensionalizeGravity(io.GetInputFileParser().ReadGravity())),
    cfl_number_(io.GetInputFileParser().ReadCFLnumber()),
    user_input_fluids_ (io.GetInputFileParser().ReadParameterOfAllFluids()),
    treated_fluids_(MapUserInputToMaterialType(user_input_fluids_)),
    fixed_boundary_values_(ComputeFixedValueBoundaryConservatives(io.GetInputFileParser())),
    logger_(logger)
{

    //Sanity Checks - Do not claim completeness.
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank == 0) {
        if(end_time_ <= 0) {throw std::invalid_argument("End time must be greater zero");}
        if(maximum_level_ > 13) {throw std::invalid_argument("Maximum level must not exceed 13");}
        if(x_number_of_level_zero_blocks_ + y_number_of_level_zero_blocks_ + z_number_of_level_zero_blocks_ == 0) {
            throw std::invalid_argument("At least one Block needs to be present on level zero");
        }
        if(block_size_on_level_zero_ <= 0) {throw std::invalid_argument("Block sizes must be greater zero");}
        if(epsilon_reference_ <= 0) {throw std::invalid_argument("Epsilon_Reference must be greater zero");}
        if(cfl_number_ < 0) {throw std::invalid_argument("CFL-Number must be larger zero");}
        if(x_number_of_level_zero_blocks_ > 128 || y_number_of_level_zero_blocks_ > 128 || z_number_of_level_zero_blocks_ > 128) {
            throw std::invalid_argument("Block number on level zero 'Block Ratio' must not exceed 128");
        }
    }
    //Logging the input parameters - Does not claim completeness.
    logger_.LogMessage("Dimensions = " + std::to_string(static_cast<int>(CC::DIM())),true,true);
    #ifdef HILBERT
    logger_.LogMessage("Hilbert Curve Load Balancing",true,true);
    #endif
    logger_.LogMessage("Starting Time t = 0.0",true,true);
    logger_.LogMessage("End Time: t = " + LogWriter::ConvertDoubleFormat(end_time_),true,true);
    logger_.LogMessage("Time naming factor: " + LogWriter::ConvertDoubleFormat(time_naming_factor_),true,true);
    logger_.LogMessage("Maximum Level: " + std::to_string(maximum_level_),true,true);
    logger_.LogMessage("CFL Number: " + std::to_string(cfl_number_),true,true);

    //NH 2017-02-09: Dirty Hack with casts ...
    std::string boundaries_string = "Boundaries: " + std::to_string(static_cast<int>(boundary_conditions_[0])) + ","
                                                   + std::to_string(static_cast<int>(boundary_conditions_[1]));
    std::string domain_size_string = "Domain Size: " + std::to_string(x_number_of_level_zero_blocks_ * block_size_on_level_zero_);
    std::string resolution_on_level_zero_string = "Resolution on Level Zero: " + std::to_string(x_number_of_level_zero_blocks_ * CC::ICX());
    std::string resolution_on_level_max_string  = "Resolution on Maximum Level: " + std::to_string(x_number_of_level_zero_blocks_ * CC::ICX() * (1 << maximum_level_));
    if(CC::DIM() != Dimension::One){
        boundaries_string.append("," + std::to_string(static_cast<int>(boundary_conditions_[2])) +
                                 "," + std::to_string(static_cast<int>(boundary_conditions_[3]))  );
        domain_size_string.append("x" + std::to_string(y_number_of_level_zero_blocks_ * block_size_on_level_zero_) );
        resolution_on_level_zero_string.append("x" + std::to_string(y_number_of_level_zero_blocks_ * CC::ICY()));
        resolution_on_level_max_string.append("x" + std::to_string(y_number_of_level_zero_blocks_ * CC::ICY() * (1 << maximum_level_)));
    }

    if(CC::DIM() == Dimension::Three){
        boundaries_string.append( "," + std::to_string(static_cast<int>(boundary_conditions_[4])) +
                                  "," + std::to_string(static_cast<int>(boundary_conditions_[5]))  );
        domain_size_string.append("x" + std::to_string(z_number_of_level_zero_blocks_ * block_size_on_level_zero_) );
        resolution_on_level_zero_string.append("x" + std::to_string(z_number_of_level_zero_blocks_ * CC::ICZ()));
        resolution_on_level_max_string.append("x" + std::to_string(z_number_of_level_zero_blocks_ * CC::ICZ() * (1 << maximum_level_)));
    }
    resolution_on_level_zero_string.append(" Internal Cells");
    resolution_on_level_max_string.append(" Internal Cells");

    logger_.LogMessage(boundaries_string, true,true);
    logger_.LogMessage(domain_size_string,true,true);
    logger_.LogMessage(resolution_on_level_zero_string);
    logger_.LogMessage("Cell Size on Level Zero: " + std::to_string(block_size_on_level_zero_/(CC::ICX())));
    logger_.LogMessage(resolution_on_level_max_string);
    //choose direction with maximum cell number
    logger_.LogMessage("Cell Size on Maximum Level: " + std::to_string(block_size_on_level_zero_/(CC::ICX() * (1 << maximum_level_))));
    //NH 2017-02-09 Dirty Hack with casting ...
    for(const auto& fluid : FluidDataForMaterialManager()) {
        logger_.LogMessage("Fluid: " + std::to_string(static_cast<int>(std::get<0>(fluid))) + " -> " + std::to_string(static_cast<int>(std::get<1>(fluid))),true,true);
        logger_.LogMessage("Gamma:" + std::to_string(std::get<2>(fluid)[0]),true,true);
        //...
    }
}

/**
 * @brief Gives the output times to be used during the simulation in treated form.
 * @param parser Inputfile parser to get the user inputs.
 * @return List of time instances at which output should be written.
 */
std::vector<double> SimulationSetup::DetermineOutputTimes(const InputFileParser& parser) const {

  std::string output_type = parser.ReadOutputTimesType();
  std::vector<double> output_times;

  double start_time = parser.ReadStartTime();
  double end_time   = parser.ReadEndTime();

  if(!(output_type.compare("Interval"))) {
    double period = parser.ReadOutputPeriod();
    if(period <= 0.0) {
        throw std::invalid_argument("Output Period must not be smaller equal zero");
    }
    double time = start_time;
    while(time < end_time) {
        time += period;
        output_times.push_back(unit_handler_.NonDimensionalizeTime(time));
    }
    if(output_times.back() >= end_time) {
        output_times.pop_back(); // We erase the last element
    }
  } else if (!output_type.compare("Timestamps")) {
        output_times = parser.ReadTimestamps();
        for(auto& time : output_times) {
            time = unit_handler_.NonDimensionalizeTime(time);
        }
        while(output_times.front() <= unit_handler_.NonDimensionalizeTime(start_time)) {
            output_times.erase(output_times.begin());
        }
        while(output_times.back() >= unit_handler_.NonDimensionalizeTime(end_time)) {
            output_times.pop_back();
        }
  } else {
    throw std::invalid_argument("Outputtype does not exist");
  }

  output_times.push_back(unit_handler_.NonDimensionalizeTime(end_time));

  // final checks
  double last_element = 0.0;
  for(const auto& time : output_times) {
    if(last_element > time) {
        throw std::invalid_argument("Timestamps must be ascending");
    }
    if(time < 0.0) {
        throw std::invalid_argument("Timestamps must be larger zero");
    }
    if(time > end_time) {
        throw std::invalid_argument("Timestamp bigger than end_time");
    }
    last_element = time;
  }

  return output_times;
}


/**
 * @brief Reads the prime states given in the inputfile for the fixed value BCs
 * and converts them into conservative states.
 * @param parser Inputfile parser to get the user inputs.
 * @return Array with conservative states in potentially all external boundary conditions.
 */
std::array<std::array<double, CC::NoEq()>, 6> SimulationSetup::ComputeFixedValueBoundaryConservatives(const InputFileParser& parser) const {

  //initialize output and temporary buffers; uses maximum number of boundary conditions (6) in 3D as this
  //makes functions below for input-file parser and conversion to primes simple
  std::array<std::array<double, CC::NoEq()>, 6> fixed_boundary_values;
  std::array<double, CC::NoEq()> temp;
  for (unsigned int i=0; i<CC::NoEq();++i) {temp[i] = 0.0;}

  //initialize fixed values at boundary with 0
  for (unsigned int i= 0; i<6; ++i){
      for (unsigned int j= 0; j<CC::NoEq(); ++j){
          fixed_boundary_values[i][j] = 0.0;
      }
  }

  //read in fixed value boundary conditions (prime states) from inputfile and convert into conservatives
  if (boundary_conditions_[0] == BoundaryType::eFixedValue){
      temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::eEast);
      fixed_boundary_values[0] = ConvertPrimesToConservatives(temp);
  }

  if (boundary_conditions_[1] == BoundaryType::eFixedValue) {
      temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::eWest);
      fixed_boundary_values[1] = ConvertPrimesToConservatives(temp);
  }

  if (CC::DIM() != Dimension::One){
      if (boundary_conditions_[2] == BoundaryType::eFixedValue){
          temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::eNorth);
          fixed_boundary_values[2] = ConvertPrimesToConservatives(temp);
      }

      if (boundary_conditions_[3] == BoundaryType::eFixedValue){
          temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::eSouth);
          fixed_boundary_values[3] = ConvertPrimesToConservatives(temp);
      }
  }


  if (CC::DIM() == Dimension::Three){
      if (boundary_conditions_[4] == BoundaryType::eFixedValue){
          temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::eTop);
          fixed_boundary_values[4] = ConvertPrimesToConservatives(temp);
      }

      if (boundary_conditions_[5] == BoundaryType::eFixedValue){
          temp = parser.ReadFixedValueBoundaryCondition(BoundaryLocation::eBottom);
          fixed_boundary_values[5] = ConvertPrimesToConservatives(temp);
      }
  }


  return fixed_boundary_values;
}


/**
 * @brief Gets prime states from inputfile and converts them into conservatives
 * required for the BCs.
 * @param values Array of fixed value BC prime states.
 * @return conservatives.
 */
std::array<double, CC::NoEq()> SimulationSetup::ConvertPrimesToConservatives(const std::array<double, CC::NoEq()> values) const {

    std::array<double, 5> conservatives = {0.0, 0.0, 0.0, 0.0, 0.0};

    //temporary material manager here necessary to compute fixed value conservatives. Only use here, generally not a good idea!!
    const MaterialManager temporary_material_manager = MaterialManager(FluidDataForMaterialManager());

    //set conservatives from inputfile prime states
    conservatives[0] = unit_handler_.NonDimensionalizeDensity(values[0]);
    conservatives[2] = unit_handler_.NonDimensionalizeVelocity(values[2])*conservatives[0];
    if (CC::DIM() != Dimension::One)
        conservatives[3] = unit_handler_.NonDimensionalizeVelocity(values[3])*conservatives[0];
    if (CC::DIM() == Dimension::Three)
        conservatives[4] = unit_handler_.NonDimensionalizeVelocity(values[4])*conservatives[0];
    conservatives[1] = temporary_material_manager.GetEnergy(treated_fluids_[0], conservatives[0],
                                                                                conservatives[2],
                                                                                conservatives[3],
                                                                                conservatives[4],
                                                                                unit_handler_.NonDimensionalizePressure(values[1]));

    //cut array to fit the return parameter
    std::array<double, CC::NoEq()> return_conservatives;
    for (unsigned int i=0; i<CC::NoEq();++i) {
        return_conservatives[i] = conservatives[i];
    }

    return return_conservatives;
}




/**
 * @brief Gives the type of the external boundary condition at the specified location.
 * @param location Direction of the edge of the domain.
 * @return The boundary identifer at the edge of the domain.
 */
BoundaryType SimulationSetup::ExternalBoundary(const BoundaryLocation location) const {

  switch (location) {
    case BoundaryLocation::eEast:
        return boundary_conditions_[0];
    break;
    case BoundaryLocation::eWest:
        return boundary_conditions_[1];
    break;
    case BoundaryLocation::eNorth:
        return boundary_conditions_[2];
    break;
    case BoundaryLocation::eSouth:
        return boundary_conditions_[3];
    break;
    case BoundaryLocation::eTop:
        return boundary_conditions_[4];
    break;
    case BoundaryLocation::eBottom:
        return boundary_conditions_[5];
    break;
    default:
        throw std::invalid_argument("Broke Boundary Conditions in Input");
    break;
    }
}

/**
 * @brief Gives the type of the external boundary condition at the specified location.
 * @param location Direction of the edge of the domain.
 * @return The boundary identifer at the edge of the domain.
 */
std::array<double,CC::NoEq()> SimulationSetup::FixedValueBoundary(const BoundaryLocation location) const {

  switch (location) {
    case BoundaryLocation::eEast:
        return fixed_boundary_values_[0];
    break;
    case BoundaryLocation::eWest:
        return fixed_boundary_values_[1];
    break;
    case BoundaryLocation::eNorth:
        return fixed_boundary_values_[2];
    break;
    case BoundaryLocation::eSouth:
        return fixed_boundary_values_[3];
    break;
    case BoundaryLocation::eTop:
        return fixed_boundary_values_[4];
    break;
    case BoundaryLocation::eBottom:
        return fixed_boundary_values_[5];
    break;
    default:
        throw std::invalid_argument("Fixed value location not found");
    break;
    }
}


/**
 * @brief Converts the fluid input from the user into input used by the MaterialManager class.
 * @return Tuple holding the Equation of State identifer, the unique material identifier and the fluid data.
 */
std::vector<std::tuple<MaterialName, MaterialName, std::vector<double>>> SimulationSetup::FluidDataForMaterialManager() const {

  std::vector<std::tuple<MaterialName, MaterialName, std::vector<double>>> result;
  if(user_input_fluids_.size() != treated_fluids_.size()) {
    throw std::invalid_argument("Fluid type vectors do not have the same size");
  }

  for(unsigned int i = 0; i < user_input_fluids_.size(); ++i) {
    result.emplace_back(std::make_tuple(user_input_fluids_[i].first, treated_fluids_[i],user_input_fluids_[i].second));
  }

  return result;
}

/**
 * @brief Gives a list of all levels which may be present in this simulation run.
 * @return Vector holding the levels in ascending order.
 */
std::vector<unsigned int> SimulationSetup::AllLevels() const {

  std::vector<unsigned int> all_levels(maximum_level_+1); //Level zero need to be counted as well
  std::iota(all_levels.begin(),all_levels.end(),0);
  return all_levels;
}

/**
 * @brief Gives the initial density at the provided location for the given material.
 * @param x_coordinates The X-Coordinates of the cell centers.
 * @param y_coordinates The Y-Coordinates of the cell centers.
 * @param z_coordinates The Z-Coordinates of the cell centers.
 * @param material The material in the block to be filled with the returned data.
 * @param initial_values Reference to array holding the resulting density. Indirect return value.
 */
void SimulationSetup::GetInitialDensity(const std::array<double,CC::ICX()>& x_coordinates,
                                        const std::array<double,CC::ICY()>& y_coordinates,
                                        const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material,
                                        double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const {

  initial_condition_.GetInitialDensity(x_coordinates,y_coordinates,z_coordinates,material,initial_values);
  for(unsigned int i = 0; i < CC::ICX(); ++i) {
    for(unsigned int j = 0; j < CC::ICY(); ++j) {
        for(unsigned int k = 0; k < CC::ICZ(); ++k) {
            initial_values[i][j][k] = unit_handler_.NonDimensionalizeDensity(initial_values[i][j][k]);
        }
    }
  }
}

/**
 * @brief Gives the initial X-Momentum at the provided location for the given material.
 * @param x_coordinates The X-Coordinates of the cell centers.
 * @param y_coordinates The Y-Coordinates of the cell centers.
 * @param z_coordinates The Z-Coordinates of the cell centers.
 * @param material The material in the block to be filled with the returned data.
 * @param rho The density at locations for which the meomenutm is to be determined.
 * @param initial_values Reference to array holding the resulting x-momentum. Indirect return value.
 */
void SimulationSetup::GetInitialMomentumX(const std::array<double,CC::ICX()>& x_coordinates,
                                          const std::array<double,CC::ICY()>& y_coordinates,
                                          const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material,
                                          double (&rho)[CC::ICX()][CC::ICY()][CC::ICZ()], double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const{

  initial_condition_.GetInitialVelocityX(x_coordinates,y_coordinates,z_coordinates,material,initial_values);
  for(unsigned int i = 0; i < CC::ICX(); ++i) {
    for(unsigned int j = 0; j < CC::ICY(); ++j) {
        for(unsigned int k = 0; k < CC::ICZ(); ++k) {
            initial_values[i][j][k] = rho[i][j][k] * unit_handler_.NonDimensionalizeVelocity(initial_values[i][j][k]);
        }
    }
  }
}

/**
 * @brief Gives the initial Y-Momentum at the provided location for the given material.
 * @param x_coordinates The X-Coordinates of the cell centers.
 * @param y_coordinates The Y-Coordinates of the cell centers.
 * @param z_coordinates The Z-Coordinates of the cell centers.
 * @param material The material in the block to be filled with the returned data.
 * @param rho The density at locations for which the meomenutm is to be determined.
 * @param initial_values Reference to array holding the resulting y-momentum. Indirect return value.
 */
void SimulationSetup::GetInitialMomentumY(const std::array<double,CC::ICX()>& x_coordinates,
                                          const std::array<double,CC::ICY()>& y_coordinates,
                                          const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material,
                                          double (&rho)[CC::ICX()][CC::ICY()][CC::ICZ()], double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const{

  initial_condition_.GetInitialVelocityY(x_coordinates,y_coordinates,z_coordinates,material,initial_values);
  for(unsigned int i = 0; i < CC::ICX(); ++i) {
    for(unsigned int j = 0; j < CC::ICY(); ++j) {
        for(unsigned int k = 0; k < CC::ICZ(); ++k) {
            initial_values[i][j][k] = rho[i][j][k] * unit_handler_.NonDimensionalizeVelocity(initial_values[i][j][k]);
        }
    }
  }

}

/**
 * @brief Gives the initial Z-Momentum at the provided location for the given material.
 * @param x_coordinates The X-Coordinates of the cell centers.
 * @param y_coordinates The Y-Coordinates of the cell centers.
 * @param z_coordinates The Z-Coordinates of the cell centers.
 * @param material The material in the block to be filled with the returned data.
 * @param rho The density at locations for which the meomenutm is to be determined.
 * @param initial_values Reference to array holding the resulting z-momentum. Indirect return value.
 */
void SimulationSetup::GetInitialMomentumZ(const std::array<double,CC::ICX()>& x_coordinates,
                                          const std::array<double,CC::ICY()>& y_coordinates,
                                          const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material,
                                          double (&rho)[CC::ICX()][CC::ICY()][CC::ICZ()],double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const{

  initial_condition_.GetInitialVelocityZ(x_coordinates,y_coordinates,z_coordinates,material,initial_values);
  for(unsigned int i = 0; i < CC::ICX(); ++i) {
    for(unsigned int j = 0; j < CC::ICY(); ++j) {
        for(unsigned int k = 0; k < CC::ICZ(); ++k) {
            initial_values[i][j][k] = rho[i][j][k] * unit_handler_.NonDimensionalizeVelocity(initial_values[i][j][k]);
        }
    }
  }

}

/**
 * @brief Gives the initial pressure at the provided location for the given material.
 * @param x_coordinates The X-Coordinates of the cell centers.
 * @param y_coordinates The Y-Coordinates of the cell centers.
 * @param z_coordinates The Z-Coordinates of the cell centers.
 * @param material The material in the block to be filled with the returned data.
 * @param initial_values Reference to array holding the resulting pressure. Indirect return value.
 */
void SimulationSetup::GetInitialPressure(const std::array<double,CC::ICX()>& x_coordinates,
                                         const std::array<double,CC::ICY()>& y_coordinates,
                                         const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material,
                                         double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const {

  initial_condition_.GetInitialPressure(x_coordinates,y_coordinates,z_coordinates,material,initial_values);
  for(unsigned int i = 0; i < CC::ICX(); ++i) {
    for(unsigned int j = 0; j < CC::ICY(); ++j) {
        for(unsigned int k = 0; k < CC::ICZ(); ++k) {
            initial_values[i][j][k] = unit_handler_.NonDimensionalizePressure(initial_values[i][j][k]);
        }
    }
  }

}

/**
 * @brief Gives the initial fluids present in the cube spanning over the provided coordinates.
 * @param x_coordinates The X-Coordinates of the cell centers.
 * @param y_coordinates The Y-Coordinates of the cell centers.
 * @param z_coordinates The Z-Coordinates of the cell centers.
 * @return An identifier for every material present in the spanned cube.
 */
std::vector<MaterialName> SimulationSetup::GetInitialFluids(const std::array<double,CC::ICX()>& x_coordinates,
                                                            const std::array<double,CC::ICY()>& y_coordinates,
                                                            const std::array<double,CC::ICZ()>& z_coordinates) const {

  // avoid compiler warnings:
  (void)x_coordinates;
  (void)y_coordinates;
  (void)z_coordinates;
  // just return the first fluid as multi-phase setups are not supported yet
  return { treated_fluids_[0] };
}

/**
 * @brief Translates the non-dimensional density values to a value with units.
 * @param rho Non-dimensional Density.
 * @return Density in unit representation.
 */
double SimulationSetup::DimensionalizeRho(const double rho) const {
    return unit_handler_.DimensionalizeDensity(rho);
}

/**
 * @brief Translates the non-dimensional momentums value to a value with units.
 * @param rho Unitless Momentum.
 * @return Momentum in unit representation.
 */
double SimulationSetup::DimensionalizeMomentum(const double momentum) const {
    return unit_handler_.DimensionalizeMomentum(momentum);
}

/**
 * @brief Translates the non-dimensional Energy value to a value with units.
 * @param rho Unitless energy.
 * @return Energy in unit representation.
 */
double SimulationSetup::DimensionalizeEnergy(const double energy) const {
    return unit_handler_.DimensionalizeEnergy(energy);
}

/**
 * @brief Translates the non-dimensional length/size value to a value with units.
 * @param rho Unitless length.
 * @return Length/size in unit representation.
 */
double SimulationSetup::DimensionalizeLength(const double length) const {
    return unit_handler_.DimensionalizeLength(length);
}

/**
 * @brief Translates the non-dimensional time value to a value with unitz.
 * @param time Unitless time.
 * @return Time in unit representation.
 */
double SimulationSetup::DimensionalizeTime(const double time) const {
    return unit_handler_.DimensionalizeTime(time);
}

/**
 * @brief Gives the smallest cell size possible in the current simulation (number of internal cells and maximum level).
 * @return .
 */
double SimulationSetup::SmallestPossibleCellSize() const {
  double level_zero_cell_size = block_size_on_level_zero_ / double(CC::ICX()); //ICX is always filled.
  double factor =  (1 << maximum_level_); // bit shift is of type "(unsigned?) int" then converison to double happens.
  return level_zero_cell_size/factor;
}

/**
 * @brief Converts the User Input Material Type, e.g. generic equation of state to a unique material
 *        identifier. Does not convert final material classes, e.g. "Water".
 * @param input The material data as retrieved from the inputfile.
 * @return Vector holding one material identifier for each input.
 */
std::vector<MaterialName> SimulationSetup::MapUserInputToMaterialType(const std::vector<std::pair<MaterialName, std::vector<double>>>& input) {

  std::vector<MaterialName> result;
  MaterialName stiffened_gas_counter = MaterialName::eUserStiffenedGasOne;
  MaterialName waterlike_fluid_counter = MaterialName::eUserWaterlikeFluidOne;
  MaterialName temp;

  for(const auto& fluid : input) {
    switch(fluid.first) {
        case MaterialName::eStiffenedGas : {
            switch (stiffened_gas_counter) {
                case MaterialName::eUserStiffenedGasOne: {
                    temp = stiffened_gas_counter;
                    stiffened_gas_counter = MaterialName::eUserStiffenedGasTwo;
                }
                break;
                case MaterialName::eUserStiffenedGasTwo: {
                    temp = stiffened_gas_counter;
                    stiffened_gas_counter = MaterialName::eUserStiffenedGasThree;
                }
                break;
                case MaterialName::eUserStiffenedGasThree: {
                    temp = stiffened_gas_counter;
                    stiffened_gas_counter = MaterialName::eStiffenedGasOutOfBounds;
                }
                break;
                default: {
                    throw std::invalid_argument("Out of Stiffened Gasses - Please use less User Defined Ones");
                }
                break;
            }
        }
        break;
        case MaterialName::eWaterlikeFluid : {
            switch (waterlike_fluid_counter){
                case MaterialName::eUserWaterlikeFluidOne: {
                temp = waterlike_fluid_counter;
                waterlike_fluid_counter = MaterialName::eUserWaterlikeFluidTwo;
            }
            break;
            case MaterialName::eUserWaterlikeFluidTwo: {
                temp = waterlike_fluid_counter;
              waterlike_fluid_counter = MaterialName::eUserWaterlikeFluidThree;
            }
            break;
            case MaterialName::eUserWaterlikeFluidThree: {
               temp = waterlike_fluid_counter;
              waterlike_fluid_counter = MaterialName::eWaterlikeFluidOutOfBounds;
            }
            break;
            default: {
                throw std::invalid_argument("Out of Waterlike Fluids - Please use less User Defined Ones");
            }
            break;
            }
        }
        break;
        default: {
            throw std::logic_error("The provided Material Type Input does not exist");
        break;
        }

    }
    result.emplace_back(temp);
  }

  return result;
}
