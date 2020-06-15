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
#include <catch.hpp>
#include "materials/material.h"
#include "materials/material_definitions.h"

#include "materials/equations_of_state/stiffened_gas.h"
#include "materials/material_property_models/shear_viscosity_models/constant_shear_viscosity_model.h"
#include "materials/material_property_models/thermal_conductivity_models/constant_thermal_conductivity_model.h"

SCENARIO( "Material is constructed and its properties are queried", "[1rank]" ) {

   // Initialize the unit handler class
   UnitHandler const unit_handler( 1.0, 1.0, 1.0, 1.0 );

   GIVEN( "A material with non-zero properties" ) {

      WHEN( "All parameters are set to non-zero values" ) {

         // Define material properties and initialize material
         double const shear_viscosity = 1.0;
         double const bulk_viscosity = 2.0;
         double const specific_heat_capacity = 3.0;
         double const thermal_heat_conductivity = 4.0;

         // Instantiate a random eos and parameter model to check if the pointer movement works properly
         std::unordered_map<std::string, double> map_data = { { "gamma", 0.0 }, { "backgroundPressure", 0.0 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGas const>( map_data, unit_handler ) );
         map_data = { { "muConstant", 0.0 } };
         std::unique_ptr<MaterialParameterModel const> shear_model( std::make_unique<ConstantShearViscosityModel const>( map_data, unit_handler ) );
         map_data = { { "lambdaConstant", 0.0 } };
         std::unique_ptr<MaterialParameterModel const> conductivity_model( std::make_unique<ConstantThermalConductivityModel const>( map_data, unit_handler ) );

         // Instantiate material
         Material const material( std::move( equation_of_state ), bulk_viscosity, shear_viscosity, thermal_heat_conductivity, specific_heat_capacity,
                                  std::move( shear_model ), std::move( conductivity_model ), unit_handler );

         THEN( "The viscosity is" ) {
            std::vector<double> expected( { 1.0, 2.0 } );
            REQUIRE( material.GetShearAndBulkViscosity() == expected );
            REQUIRE( material.GetShearViscosity() == expected.front() );
            REQUIRE( material.GetBulkViscosity() == expected.back() );
         }

         THEN( "The specific heat is" ) {
            REQUIRE( material.GetSpecificHeatCapacity() == 3.0 );
         }

         THEN( "The thermal conductivity is" ) {
            REQUIRE( material.GetThermalConductivity() == 4.0 );
         }

         THEN( "The original equation of state is" ) {
            REQUIRE( equation_of_state == nullptr );
         }

         THEN( "The material equation of state is" ) {
            REQUIRE( &material.GetEquationOfState() );
         }

         THEN( "The original shear viscosity model is" ) {
            REQUIRE( shear_model == nullptr );
         }

         THEN( "The material shear viscosity model is" ) {
            REQUIRE( &material.GetShearViscosityModel() );
         }

         THEN( "The original thermal conductivity model is" ) {
            REQUIRE( conductivity_model == nullptr );
         }

         THEN( "The material thermal conductivity model is" ) {
            REQUIRE( &material.GetThermalConductivityModel() );
         }
      }
   }
}
