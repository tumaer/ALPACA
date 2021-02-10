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
#include "materials/equations_of_state/waterlike_fluid.h"
#include "utilities/mathematical_functions.h"
#include "utilities/string_operations.h"
#include "utilities/helper_functions.h"
#include <cmath>

/**
 * @brief Constructs a waterlike equation of state with eos parameters given as input.
 * @param dimensional_eos_data Map containing all data for the equation of state.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 *
 * @note During the constructing a check is done if the required parameter exists. If not an error is thrown.
 *       Furthermore, dimensionalization of each value is done.
 */
WaterlikeFluid::WaterlikeFluid( std::unordered_map<std::string, double> const& dimensional_eos_data, UnitHandler const& unit_handler ) : gamma_( GetCheckedParameter<double>( dimensional_eos_data, "gamma", "WaterlikeFluid" ) ),
                                                                                                                                         A_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_eos_data, "A", "WaterlikeFluid" ), UnitType::Pressure ) ),
                                                                                                                                         B_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_eos_data, "B", "WaterlikeFluid" ), UnitType::Pressure ) ),
                                                                                                                                         rho0_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_eos_data, "rho0", "WaterlikeFluid" ), UnitType::Density ) ) {
   /* Empty besides initializer list*/
}

/**
 * @brief Computes pressure from inputs as A - B + B * ( rho / rho0 )^gamma.
 * @param mass The mass used for the computation.
 * @return Pressure according to Tait's equation of state.
 */
double WaterlikeFluid::ComputePressure( double const mass, double const, double const, double const, double const ) const {
   return A_ - B_ + B_ * std::pow( mass / rho0_, gamma_ );
}

/**
 * @brief Gives the enthalpy for the given inputs. ( Not available for classic Tait ).
 * @return Zero. This is according to Tait's equation of state correct.
 */
double WaterlikeFluid::ComputeEnthalpy( double const, double const, double const, double const, double const ) const {
   return 0.0;
}

/**
 * @brief Computes energy according to Taits equation.
 * @param density The density used for the computation.
 * @param velocity_x The velocity in x-direction used for the computation.
 * @param velocity_y The velocity in y-direction used for the computation.
 * @param velocity_z The velocity in z-direction used for the computation.
 * @param pressure The pressure used for the computation.
 * @return Energy according to given inputs.
 */
double WaterlikeFluid::ComputeEnergy( double const density, double const velocity_x, double const velocity_y, double const velocity_z, double const pressure ) const {
   return ( 1.0 / ( 1.0 - gamma_ ) ) * ( pressure + B_ - A_ ) + B_ - A_ + ( 0.5 * DimensionAwareConsistencyManagedSum( velocity_x * velocity_x, velocity_y * velocity_y, velocity_z * velocity_z ) * density );
}

/**
 * @brief Computes Gruneisen coefficient as ( gamma-1 ) for stiffened-gas equation of state.
 * @return Gruneisen coefficient according to Tait's equation of state.
 */
double WaterlikeFluid::GetGruneisen() const {
   return 0.0;
}

/**
 * @brief Computes psi from inputs as gamma * ( p + B - A ) / rho.
 * @param pressure The pressure used for the computation.
 * @param one_density The density used for the computation.
 * @return Psi according to Tait's equation of state.
 */
double WaterlikeFluid::ComputePsi( double const pressure, double const one_density ) const {
   return gamma_ * ( pressure + B_ - A_ ) * one_density;
}

/**
 * @brief Computes speed of sound from inputs as sqrt( gamma * ( p + B - A ) / rho ).
 * @param density The density used for the computation.
 * @param pressure The pressure used for the computation.
 * @return Speed of sound according to Tait's equation of state.
 */
double WaterlikeFluid::ComputeSpeedOfSound( double const density, double const pressure ) const {
   return std::sqrt( gamma_ * ( pressure + B_ - A_ ) / density );
}

/**
 * @brief Provides logging information of the equation of state.
 * @param indent Number of white spaces used at the beginning of each line for the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string WaterlikeFluid::GetLogData( unsigned int const indent, UnitHandler const& unit_handler ) const {
   // string that is returned
   std::string log_string;
   // Name of the equation of state
   log_string += StringOperations::Indent( indent ) + "Type                 : Water-like fluid\n";
   // Parameters with small indentation
   log_string += StringOperations::Indent( indent ) + "Gruneisen coefficient: " + StringOperations::ToScientificNotationString( GetGruneisen(), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "Gamma                : " + StringOperations::ToScientificNotationString( gamma_, 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "A                    : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( A_, UnitType::Pressure ), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "B                    : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( B_, UnitType::Pressure ), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "rho0                 : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( rho0_, UnitType::Density ), 9 ) + "\n";

   return log_string;
}
