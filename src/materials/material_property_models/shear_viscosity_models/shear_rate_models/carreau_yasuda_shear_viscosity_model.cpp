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
#include "materials/material_property_models/shear_viscosity_models/shear_rate_models/carreau_yasuda_shear_viscosity_model.h"

#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"

/**
 * @brief Constructor to instantiate the Carreau-Yasuda shear viscosity model. The model is implemented based on \cite{Bird1987}.
 *                                                                               mu0/mu-inf = viscosity at zero/infinite shear rates.
 *        mu = mu_inf + (mu0  - mu_inf) * (1 + (gamma_dot / k)^a)^((n-1)/a)          k      = the shear rate where the viscosity is constant (depending on n and a).
 *                                                                                gamma_dot = shear rate.
 *                                                                                   n      = power-law exponent (slope in power law region).
 *                                                                                   a      = transition factor between zero-shear-rate and power-law region.
 *                                                                                            A transition factor of 2 recovers the original Carreau model.
 *
 * @param dimensional_parameter_map Map, where all parameters are stored.
 * @param unit_handler Instance to provide dimensionalization of values.
 *
 * @note Runtime a check is done that all required parameters are present. Furthermore, pre-calculations are done for simpler access during compute call.
 */
CarreauYasudaShearViscosityModel::CarreauYasudaShearViscosityModel( std::unordered_map<std::string, double> const& dimensional_parameter_map,
                                                                    UnitHandler const& unit_handler ) :// Start initializer list
                                                                                                        ShearRateMaterialParameterModel<CarreauYasudaShearViscosityModel>(),
                                                                                                        mu_infinite_shear_rates_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_parameter_map, "muInfiniteShearRates", "CarreauYasudaShearViscosityModel" ), UnitType::Viscosity ) ),
                                                                                                        mu_zero_shear_rates_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_parameter_map, "muZeroShearRates", "CarreauYasudaShearViscosityModel" ), UnitType::Viscosity ) ),
                                                                                                        transition_factor_( GetCheckedParameter<double>( dimensional_parameter_map, "transitionFactor", "CarreauYasudaShearViscosityModel" ) ),
                                                                                                        power_law_exponent_( GetCheckedParameter<double>( dimensional_parameter_map, "powerLawExponent", "CarreauYasudaShearViscosityModel" ) ),
                                                                                                        shear_rate_constant_mu_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_parameter_map, "shearRateConstantMu", "CarreauYasudaShearViscosityModel" ), {}, { UnitType::Time } ) ),
                                                                                                        mu_zero_minus_infinite_( mu_zero_shear_rates_ - mu_infinite_shear_rates_ ),
                                                                                                        exponent_( ( power_law_exponent_ - 1.0 ) / transition_factor_ ),
                                                                                                        one_shear_rate_constant_mu_( 1.0 / shear_rate_constant_mu_ ) {
   /** Empty besides initializer list and friend class constructor call  */
}

/**
 * @brief Provides logging information of the given model.
 * @param indent Number of white spaces used at the beginning of each line for the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string CarreauYasudaShearViscosityModel::GetLogData( unsigned int const indent, UnitHandler const& unit_handler ) const {
   // string that is returned
   std::string log_string;
   // Add data
   log_string += StringOperations::Indent( indent ) + "Model type: Carreau-Yasuda \n";
   log_string += StringOperations::Indent( indent ) + "mu_inf    : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( mu_infinite_shear_rates_, UnitType::Viscosity ), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "mu0       : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( mu_zero_shear_rates_, UnitType::Viscosity ), 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "n         : " + StringOperations::ToScientificNotationString( power_law_exponent_, 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "a         : " + StringOperations::ToScientificNotationString( transition_factor_, 9 ) + "\n";
   log_string += StringOperations::Indent( indent ) + "constant  : " + StringOperations::ToScientificNotationString( shear_rate_constant_mu_, 9 ) + "\n";

   return log_string;
}

/**
 * @brief Computes the shear viscosity based on the model parameter and given shear rate.
 * @param shear_rate The shear rate for which the viscosity is computed.
 * @return The calculated shear viscosity.
 */
double CarreauYasudaShearViscosityModel::ComputeParameter( double const shear_rate ) const {

   double const base = 1.0 + std::pow( shear_rate * one_shear_rate_constant_mu_, transition_factor_ );
   return mu_infinite_shear_rates_ + mu_zero_minus_infinite_ * std::pow( base, exponent_ );
}
