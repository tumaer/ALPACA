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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#include "materials/equations_of_state/noble_abel_stiffened_gas.h"
#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"
#include "utilities/mathematical_functions.h"
#include <cmath>

/**
 * @brief Constructs a Nobel-Abel stiffened gas equation of state with eos parameters given as input.
 * @param dimensional_eos_data Map containing all data for the equation of state.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 *
 * @note During the constructing a check is done if the required parameter exists. If not an error is thrown.
 *       Furthermore, dimensionalization of each value is done.
 */
NobleAbelStiffenedGas::NobleAbelStiffenedGas( std::unordered_map<std::string, double> const& dimensional_eos_data, UnitHandler const& unit_handler ) : gamma_( GetCheckedParameter<double>( dimensional_eos_data, "gamma", "NobelAbelStiffendGas" ) ),
                                                                                                                                                       covolume_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_eos_data, "covolume", "NobelAbelStiffendGas" ),
                                                                                                                                                                                                       {}, { UnitType::Density } ) ),
                                                                                                                                                       pressure_constant_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_eos_data, "pressureConstant", "NobelAbelStiffendGas" ),
                                                                                                                                                                                                                UnitType::Pressure ) ),
                                                                                                                                                       energy_constant_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_eos_data, "energyConstant", "NobelAbelStiffendGas" ),
                                                                                                                                                                                                              { UnitType::Velocity, UnitType::Velocity }, {} ) ),
                                                                                                                                                       entropy_constant_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_eos_data, "entropyConstant", "NobelAbelStiffendGas" ),
                                                                                                                                                                                                               { UnitType::Velocity, UnitType::Velocity }, { UnitType::Temperature } ) ),
                                                                                                                                                       specific_heat_capacity_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_eos_data, "specificHeatCapacity", "NobelAbelStiffendGas" ),
                                                                                                                                                                                                                     { UnitType::Velocity, UnitType::Velocity }, { UnitType::Temperature } ) ) {
   /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief Computes the pressure from inputs as ( gamma - 1 ) * ( ( E  - 0.5 * rho * ||v^2|| ) / rho - e_const ) / ( 1 / rho - covolume ) - gamma * p_const.
 * @param mass The mass used for the computation.
 * @param momentum_x The momentum in x-direction used for the computation.
 * @param momentum_y The momentum in y-direction used for the computation.
 * @param momentum_z The momentum in z-direction used for the computation.
 * @param energy The energy used for the computation.
 * @return Pressure according to NASG equation of state.
 */
double NobleAbelStiffenedGas::ComputePressure( double const mass, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   double const one_mass        = 1.0 / mass;
   double const internal_energy = one_mass * ( energy - one_mass * 0.5 * DimensionAwareConsistencyManagedSum( momentum_x * momentum_x, momentum_y * momentum_y, momentum_z * momentum_z ) );

   return ( gamma_ - 1.0 ) * ( internal_energy - energy_constant_ ) / ( one_mass - covolume_ ) - gamma_ * pressure_constant_;
}

/**
 * @brief Computes enthalpy as ( E + p ) / rho.
 * @param mass The mass used for the computation.
 * @param momentum_x The momentum in x-direction used for the computation.
 * @param momentum_y The momentum in y-direction used for the computation.
 * @param momentum_z The momentum in z-direction used for the computation.
 * @param energy The energy used for the computation.
 * @return Enthalpy value.
 */
double NobleAbelStiffenedGas::ComputeEnthalpy( double const mass, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   return ( energy + ComputePressure( mass, momentum_x, momentum_y, momentum_z, energy ) ) / mass;
}

/**
 * @brief Computes Energy according to Nobel-Abel stiffend gas equation.
 * @param density The density used for the computation.
 * @param velocity_x The velocity in x-direction used for the computation.
 * @param velocity_y The velocity in y-direction used for the computation.
 * @param velocity_z The velocity in z-direction used for the computation.
 * @param pressure The pressure used for the computation.
 * @return Energy according to given inputs.
 */
double NobleAbelStiffenedGas::ComputeEnergy( double const density, double const velocity_x, double const velocity_y, double const velocity_z, double const pressure ) const {
   double const internal_energy = ( pressure + gamma_ * pressure_constant_ ) / ( density * ( gamma_ - 1.0 ) ) * ( 1.0 - covolume_ * density ) + energy_constant_;
   return density * ( internal_energy + 0.5 * DimensionAwareConsistencyManagedSum( velocity_x * velocity_x, velocity_y * velocity_y, velocity_z * velocity_z ) );
}

/**
 * @brief Computes temperature for Noble-Abel stiffened-gas EOS.
 * @param mass The mass used for the computation.
 * @param momentum_x The momentum in x-direction used for the computation.
 * @param momentum_y The momentum in y-direction used for the computation.
 * @param momentum_z The momentum in z-direction used for the computation.
 * @param energy The energy used for the computation.
 * @return Temperature according to NASG EOS.
 */
double NobleAbelStiffenedGas::ComputeTemperature( double const mass, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const {
   double const pressure = ComputePressure( mass, momentum_x, momentum_y, momentum_z, energy );
   return ( pressure + pressure_constant_ ) / ( mass * ( gamma_ - 1.0 ) * specific_heat_capacity_ ) * ( 1.0 - covolume_ * mass );
}

/**
 * @brief Computes Gruneisen coefficient as ( gamma - 1 ) / ( 1 - covolume * rho )
 * @return Gruneisen coefficient according to NASG equation of state.
 */
double NobleAbelStiffenedGas::GetGruneisen( double const density ) const {
   return ( gamma_ - 1.0 ) / ( 1.0 - covolume_ * density );
}

/**
 * @brief DO NOT CALL. Throws an error since Gruneisen coefficient for NASG is density dependent. Call NobleAbelStiffenedGas::GetGruneisen( double const density ) instead.
 */
double NobleAbelStiffenedGas::GetGruneisen() const {
   throw std::runtime_error( "NobleAbelStiffenedGas: Gruneisen parameter depends on density!" );
}

/**
 * @brief Computes psi from inputs as ( p + gamma * p_const ) / ( rho * ( 1 - covolume * rho ) ).
 * @param pressure The pressure used for the computation.
 * @param one_density ( devision by zero is already avoided before ) .
 * @return Psi according to NASG equation of state.
 */
double NobleAbelStiffenedGas::ComputePsi( double const pressure, double const one_density ) const {
   return ( pressure + gamma_ * pressure_constant_ ) * one_density * one_density / ( one_density - covolume_ );
}

/**
 * @brief Returns Gamma.
 * @return Gamma.
 */
double NobleAbelStiffenedGas::GetGamma() const {
   return gamma_;
}

/**
 * @brief Returns the pressure constant.
 * @return B.
 */
double NobleAbelStiffenedGas::GetB() const {
   return pressure_constant_;
}

/**
 * @brief Computes Speed of Sound from inputs as sqrt( gamma * ( p + p_const ) / ( rho * ( 1 - covolume * rho ) ) ).
 * @param density The density used for the computation.
 * @param pressure The pressure used for the computation.
 * @return Speed of sound according to NASG equation of state.
 */
double NobleAbelStiffenedGas::ComputeSpeedOfSound( double const density, double const pressure ) const {
   double const speed_of_sound_squared = gamma_ * ( pressure + pressure_constant_ ) / ( density * ( 1.0 - covolume_ * density ) );
   return std::sqrt( speed_of_sound_squared );
}

/**
 * @brief Provides logging information of the equation of state.
 * @param indent Number of white spaces used at the beginning of each line for the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string NobleAbelStiffenedGas::GetLogData( unsigned int const indent, UnitHandler const& unit_handler ) const {
   // string that is returned
   std::string log_string;
   // Name of the equation of state
   log_string += StringOperations::Indent( indent ) + "Type                  : Nobel-Abel stiffened gas \n";
   // Parameters with small indentation
   log_string += StringOperations::Indent( indent ) + "Gruneisen coefficient : Density dependent \n";
   log_string += StringOperations::Indent( indent ) + "Gamma                 : " + StringOperations::ToScientificNotationString( gamma_, 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "Co-volume             : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( covolume_, {}, { UnitType::Density } ), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "Pressure constant     : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( pressure_constant_, UnitType::Pressure ), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "Energy constant       : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( energy_constant_, { UnitType::Velocity, UnitType::Velocity }, {} ), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "Entropy constant      : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( entropy_constant_, { UnitType::Velocity, UnitType::Velocity }, { UnitType::Temperature } ), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "Specific heat capacity: " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( specific_heat_capacity_, { UnitType::Velocity, UnitType::Velocity }, { UnitType::Temperature } ), 9 ) + "\n";

   return log_string;
}
