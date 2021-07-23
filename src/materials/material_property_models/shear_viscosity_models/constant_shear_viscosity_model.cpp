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
#include "materials/material_property_models/shear_viscosity_models/constant_shear_viscosity_model.h"

#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"

/**
 * @brief Constructor to instantiate the constant shear viscosity model.
 * @param dimensional_parameter_map Map, where all parameters are stored.
 * @param unit_handler Instance to provide dimensionalization of values.
 *
 * @note Runtime a check is done that all required parameters are present. Furthermore, pre-calculations are done for simpler access during compute call.
 */
ConstantShearViscosityModel::ConstantShearViscosityModel( std::unordered_map<std::string, double> const& dimensional_parameter_map,
                                                          UnitHandler const& unit_handler ) :// Start initializer list
                                                                                              ConstantMaterialParameterModel<ConstantShearViscosityModel>(),
                                                                                              mu_constant_( unit_handler.NonDimensionalizeValue( GetCheckedParameter<double>( dimensional_parameter_map, "muConstant", "ConstantShearViscosityModel" ), UnitType::Viscosity ) ) {
   /** Empty besides initializer list and friend class constructor call  */
}

/**
 * @brief Provides logging information of the given model.
 * @param indent Number of white spaces used at the beginning of each line for the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string ConstantShearViscosityModel::GetLogData( unsigned int const indent, UnitHandler const& unit_handler ) const {
   // string that is returned
   std::string log_string;
   // Add data
   log_string += StringOperations::Indent( indent ) + "Model type: Constant \n";
   log_string += StringOperations::Indent( indent ) + "mu_const  : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( mu_constant_, UnitType::Viscosity ), 9 ) + "\n";

   return log_string;
}

/**
 * @brief Provides a constant shear viscosity from user input
 * @return The constant shear viscosity
 */
double ConstantShearViscosityModel::ComputeParameter() const {
   return mu_constant_;
}
