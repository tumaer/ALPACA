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

#ifndef SIMULATION_SETUP_H
#define SIMULATION_SETUP_H

#include <vector>
#include <array>
#include <tuple>

#include <sys/stat.h> // stat
#include <fstream>

#include "boundary_condition/boundary_specifications.h"
#include "enums/dimension_definition.h"
#include "materials/material_names.h"
#include "inputfile_parser.h"
#include "initial_condition.h"
#include "input_output_manager.h"
#include "unit_handler.h"
#include "output/output_types.h"
#include "log_writer.h"

/**
 * @brief The SimulationSetup class gives the interface between the user and the simulation kernel. The user input via the xml file is processed
 *        and the given data converted into quantities usable by the kernel. E.g. inputs are given in terms of primary states and are converted to
 *        conservative states.
 */
class SimulationSetup {
  const InputOutputManager& input_output_;

  //Delegating Objects
  const InitialCondition initial_condition_;
  const UnitHandler unit_handler_;

  // IO-Parameters
  const OutputType output_format_;
  const std::vector<double> output_times_;

  //time parameters
  const double start_time_;
  const double end_time_;
  const double time_naming_factor_;

  //domain parameters
  const unsigned int maximum_level_;
  const unsigned short x_number_of_level_zero_blocks_;
  const unsigned short y_number_of_level_zero_blocks_;
  const unsigned short z_number_of_level_zero_blocks_;
  const double block_size_on_level_zero_;
  const std::array<BoundaryType,6> boundary_conditions_;


  // MR parameters
  double epsilon_reference_;


  //fluid parameters
  const std::array<double, 3> gravity_;
  const double cfl_number_;
  const std::vector<std::pair<MaterialName, std::vector<double>>>  user_input_fluids_;
  const std::vector<MaterialName> treated_fluids_;

  //fixed value BC parameters - has to be defined here after the fluid parameters
  const std::array<std::array<double,CC::NoEq()>,6> fixed_boundary_values_;


  //Logger must not be const (otherwise no logbook cannot be appended
  LogWriter& logger_;

  std::vector<double> DetermineOutputTimes(const InputFileParser& parser) const;

  //functions for conversion of fixed value BC parameters from prime states into conservatives
  std::array<std::array<double, CC::NoEq()>, 6> ComputeFixedValueBoundaryConservatives(const InputFileParser& parser) const;
  std::array<double, CC::NoEq()> ConvertPrimesToConservatives(const std::array<double, CC::NoEq()> values) const;

  /**
   * @brief Computes the reference epsilon used in the multiresolution analysis, i.e. coarsening/refinement decision.
   * @param user_epsilon_reference The threshold value the user wants to ensure on the given level.
   * @param user_level_of_reference The level at which the given threshold should be ensured.
   * @return Epsilon Reference to be used in multiresolution analysis, i.e. eps_r in the equation: eps_l = 2^(D*(l-L_max)) * eps_ref.
   */
  inline double ComputeReferenceEpsilon(const double user_epsilon_reference, const unsigned int user_level_of_reference) const {
    return user_epsilon_reference *  std::pow(2, (CC::MROC()+1) * int(int(user_level_of_reference) - int(maximum_level_)) );
  }

   public:
    SimulationSetup() = delete;
    SimulationSetup(const InputOutputManager& io, LogWriter& logger);

    BoundaryType ExternalBoundary(const BoundaryLocation location) const;
    std::array<double,CC::NoEq()> FixedValueBoundary(const BoundaryLocation location) const;

    std::vector<std::tuple<MaterialName, MaterialName, std::vector<double>>> FluidDataForMaterialManager() const;

    std::vector<unsigned int> AllLevels() const;

    void GetInitialDensity(  const std::array<double,CC::ICX()>& x_coordinates, const std::array<double,CC::ICY()>& y_coordinates,const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material, double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const;
    void GetInitialPressure( const std::array<double,CC::ICX()>& x_coordinates, const std::array<double,CC::ICY()>& y_coordinates,const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material, double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const;
    void GetInitialMomentumX(const std::array<double,CC::ICX()>& x_coordinates, const std::array<double,CC::ICY()>& y_coordinates,const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material, double (&rho)[CC::ICX()][CC::ICY()][CC::ICZ()], double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const;
    void GetInitialMomentumY(const std::array<double,CC::ICX()>& x_coordinates, const std::array<double,CC::ICY()>& y_coordinates,const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material, double (&rho)[CC::ICX()][CC::ICY()][CC::ICZ()], double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const;
    void GetInitialMomentumZ(const std::array<double,CC::ICX()>& x_coordinates, const std::array<double,CC::ICY()>& y_coordinates,const std::array<double,CC::ICZ()>& z_coordinates, const MaterialName material, double (&rho)[CC::ICX()][CC::ICY()][CC::ICZ()], double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const;
    std::vector<MaterialName> GetInitialFluids(const std::array<double,CC::ICX()>& x_coordinates, const std::array<double,CC::ICY()>& y_coordinates,const std::array<double,CC::ICZ()>& z_coordinates) const;

    double DimensionalizeRho(const double rho) const;
    double DimensionalizeMomentum(const double momentum) const;
    double DimensionalizeEnergy(const double energy) const;
    double DimensionalizeLength(const double length) const;
    double DimensionalizeTime(const double time) const;

    double SmallestPossibleCellSize() const;

    /**
     * @brief Gives the time at which the simulation should start.
     * @return Dimensionless time at which the simulation should start.
     */
    inline double GetStartTime() const {return start_time_;}


    /**
     * @brief Gives the time at which the simulation should terminate.
     * @return Dimensionless time at which the simulation should end.
     */
    inline double GetEndTime() const {return end_time_;}

    /**
     * @brief Gives the dimensional multiplication factor for the output file timestamps.
     * @return Dimensional multiplication factor for output time files.
     */
    inline double GetTimeNamingFactor() const {return time_naming_factor_;}

    /**
     * @brief Gives the times at which output should be created.
     * @return List of output times.
     */
    inline std::vector<double> OutputTimes() const {return output_times_;}
    /**
     * @brief Gives the Courant–Friedrichs–Lewy number.
     * @return CFL number.
     */
    inline double GetCflNumber() const {return cfl_number_;}
    /**
     * @brief Gives the size of a block on level zero.
     * @return Block Size on level 0.
     */
    inline double LevelZeroBlockSize() const {return block_size_on_level_zero_;}
    /**
     * @brief Gives the adjusted epsilon_ref, i.e. only needs to be adjusted by the level and the dimension factor.
     * @return Epsilon_ref.
     */
    inline double EpsilonReference() const {return epsilon_reference_;}
    /**
     * @brief Gives the output format to be used.
     * @return Output format.
     */
    inline OutputType GetOutputFormat() const {return output_format_;}
    /**
     * @brief Gives the maximum level used in this simulation.
     * @return Maximum level.
     */
    inline unsigned int GetMaximumLevel() const {return maximum_level_;}
    /**
     * @brief Gives the number of blocks on level zero in X-Direction.
     * @return Number of blocks in X-Direction on Level zero.
     */
    inline unsigned short GetLevelZeroBlocksX() const {return x_number_of_level_zero_blocks_;}
    /**
     * @brief Gives the number of blocks on level zero in Y-Direction.
     * @return Number of blocks in Y-Direction on Level zero.
     */
    inline unsigned short GetLevelZeroBlocksY() const {return y_number_of_level_zero_blocks_;}
    /**
     * @brief Gives the number of blocks on level zero in Z-Direction.
     * @return Number of blocks in Z-Direction on Level zero.
     */
    inline unsigned short GetLevelZeroBlocksZ() const {return z_number_of_level_zero_blocks_;}
    /**
     * @brief Gives the gravity as array 0: X-Component, 1: Y-Component, 2: Z-Component.
     * @return Gravity.
     */
    inline std::array<double, 3> GetGravity() const {return gravity_;}

    static std::vector<MaterialName> MapUserInputToMaterialType(const std::vector<std::pair<MaterialName, std::vector<double>>>& input);
};

#endif // SIMULATION_SETUP_H
