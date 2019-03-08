/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#include "materials/equations_of_state/stiffened_gas.h"
#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"
#include "utilities/mathematical_functions.h"
#include <cmath>

#include <iostream>

/**
 * @brief Constructs a stiffened gas equation of state with eos parameters given as input.
 * @param dimensional_eos_data Map containing all data for the equation of state.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 *
 * @note During the constructing a check is done if the required parameter exists. If not an error is thrown.
 *       Furthermore, dimensionalization of each value is done.
 */
StiffenedGas::StiffenedGas( std::unordered_map<std::string, double> const& dimensional_eos_data, UnitHandler const& unit_handler ) :
   gamma_( GetCheckedParameter( dimensional_eos_data, "gamma", "StiffenedGas" ) ),
   background_pressure_( unit_handler.NonDimensionalizeValue( GetCheckedParameter( dimensional_eos_data, "B", "StiffenedGas" ), UnitType::Pressure ) ) {
   /* Empty besides initializer list*/
}

/**
 * @brief Computes Pressure from inputs as -gamma*B + ( gamma - 1 ) * ( E  - 0.5 * rho * ||v^2|| )
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Pressure according to stiffened-gas equation of state.
 */
double StiffenedGas::DoGetPressure( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   return -gamma_ * background_pressure_ + ( gamma_ - 1.0 ) * ( energy - 0.5 * DimensionAwareConsistencyManagedSum( momentum_x * momentum_x,
                                                                                                                    momentum_y * momentum_y,
                                                                                                                    momentum_z * momentum_z )
                                                                / density );
}

/**
 * @brief Computes enthalpy as ( E + p ) / rho.
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Enthalpy value.
 */
double StiffenedGas::DoGetEnthalpy( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   return ( energy + DoGetPressure( density, momentum_x, momentum_y, momentum_z, energy ) ) / density;
}

/**
 * @brief Computes Energy from inputs as ( p + gamma * B ) / ( gamma - 1 ) + 0.5 * rho * ||v^2||
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param pressure .
 * @return Energy according to stiffened-gas equation of state.
 */
double StiffenedGas::DoGetEnergy( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const pressure ) const{
   return ( pressure + gamma_ * background_pressure_ ) / ( gamma_ -1.0 )
          + ( 0.5 * DimensionAwareConsistencyManagedSum( momentum_x * momentum_x, momentum_y * momentum_y, momentum_z * momentum_z ) / density );
}

/**
 * @brief Computes Gruneisen coefficient as ( gamma-1 ) for stiffened-gas equation of state
 * @return Gruneisen coefficient .
 */
double StiffenedGas::DoGetGruneisen() const {
   return ( gamma_ - 1.0 );
}

/**
 * @brief Returns Gamma
 * @return Gamma.
 */
double StiffenedGas::DoGetGamma() const {
   return gamma_;
}

/**
 * @brief Returns B
 * @return B.
 */
double StiffenedGas::DoGetB() const {
   return background_pressure_;
}

/**
 * @brief Computes psi from inputs as ( p + gamma * B ) / rho
 * @param pressure .
 * @param one_density .
 * @return Psi according to stiffened-gas equation of state.
 */
double StiffenedGas::DoGetPsi( double const pressure, double const one_density ) const {
   return ( pressure + gamma_ * background_pressure_ ) * one_density;
}

/**
 * @brief Computes Speed of Sound from inputs as sqrt( gamma * ( p + B ) ) / rho
 * @param density .
 * @param pressure .
 * @return Speed of sound according to stiffened-gas equation of state.
 */
double StiffenedGas::DoGetSpeedOfSound( double const density, double const pressure ) const{
   return std::sqrt( gamma_*( pressure + background_pressure_ ) / density );
}

/**
 * @brief Provides logging information of the equation of state.
 * @param indent Number of white spaces used at the beginning of each line for the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string StiffenedGas::GetLogData( unsigned int const indent, UnitHandler const& unit_handler ) const {
   // string that is returned
   std::string log_string;
   // Name of the equation of state
   log_string += StringOperations::Indent( indent ) + "Type                 : Stiffened gas\n";
   // Parameters with small indentation
   log_string += StringOperations::Indent( indent ) + "Gruneisen coefficient: " + StringOperations::ToScientificNotationString(
                   DoGetGruneisen(), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "Gamma                : " + StringOperations::ToScientificNotationString(
                   gamma_, 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "Background pressure  : " + StringOperations::ToScientificNotationString(
                   unit_handler.DimensionalizeValue( background_pressure_, UnitType::Pressure ), 9 ) + "\n";
   return log_string;
}