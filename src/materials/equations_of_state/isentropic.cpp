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

#include "isentropic.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief Constructs an isentropic-eos object with material parameters given as input.
 * @param dimensional_eos_data Map containing all data for the equation of state.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @note During the constructing a check is done if the required parameter exists. If not an error is thrown.
 *       Furthermore, dimensionalization of each value is done.
 */
Isentropic::Isentropic( std::unordered_map<std::string, double> const& dimensional_eos_data, UnitHandler const& unit_handler ) : gamma_( GetCheckedParameter( dimensional_eos_data, "gamma", "Isentropic" ) ),
                                                                                                                                 A_( unit_handler.NonDimensionalizeValue( GetCheckedParameter( dimensional_eos_data, "A", "Isentropic" ), UnitType::Pressure ) ) {
   /* Empty besides initializer list*/
}

/**
 * @brief Computes pressure from inputs as A * rho^gamma.
 * @param mass The mass used in the computation.
 * @note Only the mass is needed in the computation in the isentropic case.
 */
double Isentropic::ComputePressure( double const mass, double const, double const, double const, double const ) const {
   return A_ * std::pow( mass, gamma_ );
}

/*
 * @brief Enthalpy does not matter in isentropic case.
 * @return -1.
 */
double Isentropic::ComputeEnthalpy( double const, double const, double const, double const, double const ) const {
   return -1;
}

/**
 * @brief Energy does not matter in isentropic case.
 * @return -1.
 */
double Isentropic::ComputeEnergy( double const, double const, double const, double const, double const ) const {
   return -1;
}

/**
 * @brief Gruneisen coefficent is not necessary for isentropic equation of state.
 * @return -1.
 */
double Isentropic::GetGruneisen() const {
   return -1;
}

/**
 * @brief Returns gamma.
 * @return gamma.
 */
double Isentropic::GetGamma() const {
   return gamma_;
}

/**
 * @brief Background pressure is not necessary for isentropic equation of state.
 * @return -1.
 */
double Isentropic::GetB() const {
   return -1;
}

/**
 * @brief Psi is not necessary for isentropic equation of state.
 * @param pressure .
 * @param one_density (1 / rho).
 * @return -1.
 */
double Isentropic::ComputePsi( double const, double const ) const {
   return -1;
}

/**
 * @brief Computes speed of sound from inputs as sqrt(gamma * p) / rho.
 * @param density .
 * @param pressure .
 * @return Speed of sound according to isentropic equation of state.
 */
double Isentropic::ComputeSpeedOfSound( double const density, double const pressure ) const {
   return std::sqrt( gamma_ * pressure / density );
}

/**
 * @brief Provides logging information of the equation of state.
 * @param indent Number of white spaces used at the beginning of each line for the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string Isentropic::GetLogData( unsigned int const indent, UnitHandler const& unit_handler ) const {
   // string that is returned
   std::string log_string;
   // Name of the equation of state
   log_string += StringOperations::Indent( indent ) + "Type                 : Isentropic\n";
   // Parameters with small indentation
   log_string += StringOperations::Indent( indent ) + "Gamma                : " + StringOperations::ToScientificNotationString( gamma_, 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "A                    : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( A_, UnitType::Pressure ), 9 ) + "\n";

   return log_string;
}
