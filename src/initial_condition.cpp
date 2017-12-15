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

#include "initial_condition.h"
#include "inputfile_parser.h"
#include "simulation_setup.h"

#include <stdexcept>

const std::string InitialCondition::var_density_ = "DENSITY";
const std::string InitialCondition::var_velocity_x_ = "VELOCITY_X";
const std::string InitialCondition::var_velocity_y_ = "VELOCITY_Y";
const std::string InitialCondition::var_velocity_z_ = "VELOCITY_Z";
const std::string InitialCondition::var_pressure_ = "PRESSURE";

/**
 * @brief Default constructor requiring a parse-able inputfile %Currently only XML file%
 * @param parser An user inputfile parser instance.
 */
InitialCondition::InitialCondition(const InputFileParser& parser) :
  number_of_fluids_(parser.ReadNumberOfFluids()),
  fluid_initialisation_(parser.ReadInitialConditionOfFluids()),
  fluids_(SimulationSetup::MapUserInputToMaterialType(parser.ReadParameterOfAllFluids()))
{
}

/**
 * @brief Compiles the given expression such that it is capable to give the desired fluid variables.
 * @param expression The expression in original text form.
 * @param x,y,z Reference to spatial variable.
 * @return The compiled expression ready to be evaluated.
 */
UserExpression InitialCondition::CreateFluidExpression(std::string expression, double &x, double &y, double &z) const {
  std::vector<std::tuple<std::string,double&>> variables_in;
  variables_in.push_back(std::make_tuple(std::string("x"),std::ref(x)));
  variables_in.push_back(std::make_tuple(std::string("y"),std::ref(y)));
  variables_in.push_back(std::make_tuple(std::string("z"),std::ref(z)));

  std::vector<std::string> variables_out;
  variables_out.push_back(var_density_);
  variables_out.push_back(var_velocity_x_);
  variables_out.push_back(var_velocity_y_);
  variables_out.push_back(var_velocity_z_);
  variables_out.push_back(var_pressure_);

  return UserExpression(expression,variables_in,variables_out);
}

/**
 * @brief Evaluates the user input for the initial density at the provided X-/Y- and Z-coordinates.
 * @param x_coordinates X-coordinates of the cell centers.
 * @param y_coordinates Y-coordinates of the cell centers.
 * @param z_coordinates Z-coordinates of the cell centers.
 * @param material Identifier to fill the cells with the state associated with its material.
 */
void InitialCondition::GetInitialDensity(const std::array<double, CC::ICX()>& x_coordinates,
                                         const std::array<double, CC::ICY()>& y_coordinates,
                                         const std::array<double, CC::ICZ()>& z_coordinates, MaterialName material,
                                         double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const {

  std::string Fluid_ini;

  for(int ii = 0; ii < number_of_fluids_; ++ii){
    if(fluids_[ii] == material){
        Fluid_ini = fluid_initialisation_[ii];
    }
  }

  double running_x;
  double running_y;
  double running_z;

  UserExpression density_expr = CreateFluidExpression(Fluid_ini, running_x, running_y, running_z);

  for(unsigned int i = 0; i < x_coordinates.size(); ++i) {
      running_x = x_coordinates[i];
      for(unsigned int j = 0; j < y_coordinates.size(); ++j) {
          running_y = y_coordinates[j];
          for(unsigned int k = 0; k < z_coordinates.size(); ++k) {
              running_z = z_coordinates[k];
              initial_values[i][j][k] = density_expr.GetValue(var_density_);
          }
      }
  }

}

/**
 * @brief Evaluates the user input for the initial velocity component in X-direction at the provided X-/Y- and Z-coordinates.
 * @param x_coordinates X-coordinates of the cell centers.
 * @param y_coordinates Y-coordinates of the cell centers.
 * @param z_coordinates Z-coordinates of the cell centers.
 * @param material Identifier to fill the cells with the state associated with its material.
 */
void InitialCondition::GetInitialVelocityX(const std::array<double, CC::ICX()>& x_coordinates,
                                           const std::array<double, CC::ICY()>& y_coordinates,
                                           const std::array<double, CC::ICZ()>& z_coordinates, MaterialName material,
                                           double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const {

  std::string Fluid_ini;

  for(int ii = 0; ii < number_of_fluids_; ++ii){
    if(fluids_[ii] == material){
        Fluid_ini = fluid_initialisation_[ii];
    }
  }

  double running_x;
  double running_y;
  double running_z;

  UserExpression velocity_x_expr = CreateFluidExpression(Fluid_ini, running_x, running_y, running_z);

  for(unsigned int i = 0; i < x_coordinates.size(); ++i) {
    running_x = x_coordinates[i];
    for(unsigned int j = 0; j < y_coordinates.size(); ++j) {
        running_y = y_coordinates[j];
        for(unsigned int k = 0; k < z_coordinates.size(); ++k) {
            running_z = z_coordinates[k];
            initial_values[i][j][k] = velocity_x_expr.GetValue(var_velocity_x_);
        }
    }
  }

}

/**
 * @brief Evaluates the user input for the initial velocity component in Y-direction at the provided X-/Y- and Z-coordinates.
 * @param x_coordinates X-coordinates of the cell centers.
 * @param y_coordinates Y-coordinates of the cell centers.
 * @param z_coordinates Z-coordinates of the cell centers.
 * @param material Identifier to fill the cells with the state associated with its material.
 */
void InitialCondition::GetInitialVelocityY(const std::array<double, CC::ICX()>& x_coordinates,
                                           const std::array<double, CC::ICY()>& y_coordinates,
                                           const std::array<double, CC::ICZ()>& z_coordinates, MaterialName material,
                                           double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const {

  std::string Fluid_ini;

  for(int ii = 0; ii < number_of_fluids_; ++ii){
    if(fluids_[ii] == material){
        Fluid_ini = fluid_initialisation_[ii];
    }
  }

  double running_x;
  double running_y;
  double running_z;

  UserExpression velocity_y_expr = CreateFluidExpression(Fluid_ini, running_x, running_y, running_z);

  for(unsigned int i = 0; i < x_coordinates.size(); ++i) {
    running_x = x_coordinates[i];
    for(unsigned int j = 0; j < y_coordinates.size(); ++j) {
        running_y = y_coordinates[j];
        for(unsigned int k = 0; k < z_coordinates.size(); ++k) {
            running_z = z_coordinates[k];
            initial_values[i][j][k] = velocity_y_expr.GetValue(var_velocity_y_);
        }
    }
  }
}

/**
 * @brief Evaluates the user input for the initial velocity component in Z-direction at the provided X-/Y- and Z-coordinates.
 * @param x_coordinates X-coordinates of the cell centers.
 * @param y_coordinates Y-coordinates of the cell centers.
 * @param z_coordinates Z-coordinates of the cell centers.
 * @param material Identifier to fill the cells with the state associated with its material.
 */
void InitialCondition::GetInitialVelocityZ(const std::array<double, CC::ICX()>& x_coordinates,
                                           const std::array<double, CC::ICY()>& y_coordinates,
                                           const std::array<double, CC::ICZ()>& z_coordinates, MaterialName material,
                                           double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const {

  std::string Fluid_ini;

  for(int ii = 0; ii < number_of_fluids_; ++ii){
    if(fluids_[ii] == material){
        Fluid_ini = fluid_initialisation_[ii];
    }
  }

  double running_x;
  double running_y;
  double running_z;

  UserExpression velocity_z_expr = CreateFluidExpression(Fluid_ini, running_x, running_y, running_z);

  for(unsigned int i = 0; i < x_coordinates.size(); ++i) {
    running_x = x_coordinates[i];
    for(unsigned int j = 0; j < y_coordinates.size(); ++j) {
        running_y = y_coordinates[j];
        for(unsigned int k = 0; k < z_coordinates.size(); ++k) {
            running_z = z_coordinates[k];
            initial_values[i][j][k] = velocity_z_expr.GetValue(var_velocity_z_);
        }
    }
  }
}

/**
 * @brief Evaluates the user input for the initial pressure at the provided X-/Y- and Z-coordinates.
 * @param x_coordinates X-coordinates of the cell centers.
 * @param y_coordinates Y-coordinates of the cell centers.
 * @param z_coordinates Z-coordinates of the cell centers.
 * @param material Identifier to fill the cells with the state associated with its material.
 */
void InitialCondition::GetInitialPressure(const std::array<double, CC::ICX()>& x_coordinates,
                                          const std::array<double, CC::ICY()>& y_coordinates,
                                          const std::array<double, CC::ICZ()>& z_coordinates, MaterialName material,
                                          double (&initial_values)[CC::ICX()][CC::ICY()][CC::ICZ()]) const{

  std::string Fluid_ini;

  for(int ii = 0; ii < number_of_fluids_; ++ii){
    if(fluids_[ii] == material){
        Fluid_ini = fluid_initialisation_[ii];
    }
  }

  double running_x;
  double running_y;
  double running_z;

  UserExpression pressure_expr = CreateFluidExpression(Fluid_ini, running_x, running_y, running_z);

  for(unsigned int i = 0; i < x_coordinates.size(); ++i) {
      running_x = x_coordinates[i];
      for(unsigned int j = 0; j < y_coordinates.size(); ++j) {
          running_y = y_coordinates[j];
          for(unsigned int k = 0; k < z_coordinates.size(); ++k) {
              running_z = z_coordinates[k];
              initial_values[i][j][k] = pressure_expr.GetValue(var_pressure_);
          }
      }
  }

}
