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
#include <catch2/catch.hpp>

#include <string>
#include <memory>
#include <unordered_map>

#include "input_output/input_reader/material_reader/material_reader.h"
#include "input_output/input_reader/material_reader/xml_material_reader.h"

SCENARIO( "Check that the xml material reader works properly", "[1rank]" ) {
   GIVEN( "A xml document with the valid content to read the output data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <materials>"
                                  "     <numberOfMaterials> 4 </numberOfMaterials>"
                                  "     <material1>"
                                  "       <type> SolidBounDary </type>"
                                  "       <equationOfState>"
                                  "          <type> StiffenedGas </type>"
                                  "          <type1> 0.1 </type1>"
                                  "          <param1> 0.2 </param1>"
                                  "          <param2> 1e-3 </param2>"
                                  "       </equationOfState>"
                                  "       <properties>"
                                  "          <specificHeatCapacity> 1e-4 </specificHeatCapacity>"
                                  "          <thermalConductivity> 2.0"
                                  "             <fixedValue> 1e-3 </fixedValue>"
                                  "          </thermalConductivity>"
                                  "          <bulkViscosity>"
                                  "             <fixedValue> 1e-3 </fixedValue>"
                                  "          </bulkViscosity>"
                                  "          <shearViscosity> 2.0"
                                  "            <model>"
                                  "               <name> Constant </name>"
                                  "               <parameter>"
                                  "                  <modelParam> 1e-5 </modelParam>"
                                  "               </parameter>"
                                  "            </model>"
                                  "          </shearViscosity>"
                                  "       </properties>"
                                  "     </material1>"
                                  "  </materials>"
                                  "  <materialPairings>"
                                  "     <material1_2>"
                                  "        <surfaceTensionCoefficient> 2.0"
                                  "          <model>"
                                  "             <name> Constant </name>"
                                  "             <parameter>"
                                  "                <modelParam> 1e-6 </modelParam>"
                                  "             </parameter>"
                                  "          </model>"
                                  "        </surfaceTensionCoefficient>"
                                  "     </material1_2>"
                                  "  </materialPairings>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<MaterialReader const> const reader( std::make_unique<XmlMaterialReader const>( xml_tree ) );

      WHEN( "The number of materials is read." ) {
         THEN( "The number of materials should be 4" ) {
            REQUIRE( reader->ReadNumberOfMaterials() == 4 );
         }
      }

      WHEN( "The type of materials is read." ) {
         THEN( "The type should be SolidBoundary" ) {
            REQUIRE( reader->ReadMaterialType( 1, MaterialType::Fluid ) == MaterialType::SolidBoundary );
         }
      }

      WHEN( "The equation of state data is read." ) {
         THEN( "The type should be StiffenedGas and three parameters should be in the map (type1=0.1, param1=0.2, param2=1e-3)." ) {
            std::unordered_map<std::string, double> const parameters( reader->ReadEquationOfStateData( 1 ) );
            REQUIRE( reader->ReadEquationOfStateName( 1 ) == EquationOfStateName::StiffenedGas );
            REQUIRE( parameters.at( "type1" ) == 0.1 );
            REQUIRE( parameters.at( "param1" ) == 0.2 );
            REQUIRE( parameters.at( "param2" ) == 1e-3 );
         }
      }

      WHEN( "The fixed value specific heat capacity is read." ) {
         THEN( "The fixed value should be 1e-4" ) {
            REQUIRE( reader->ReadFixedValue( { 1 }, MaterialProperty::SpecificHeatCapacity ) == 1e-4 );
         }
      }

      WHEN( "The fixed value thermal conductivity is read." ) {
         THEN( "The fixed value should be 2.0" ) {
            REQUIRE( reader->ReadFixedValue( { 1 }, MaterialProperty::ThermalConductivity ) == 2.0 );
         }
      }

      WHEN( "The fixed value bulk viscosity is read." ) {
         THEN( "The fixed value should be 1.0" ) {
            REQUIRE( reader->ReadFixedValue( { 1 }, MaterialProperty::BulkViscosity ) == 1e-3 );
         }
      }

      WHEN( "The shear viscosity model is read." ) {
         THEN( "The model name should be constant and the parameter should be modelParam=1e-5." ) {
            std::unordered_map<std::string, double> const parameters( reader->ReadModelData( { 1 }, MaterialProperty::ShearViscosity ) );
            REQUIRE( reader->ReadModelName( { 1 }, MaterialProperty::ShearViscosity ) == MaterialPropertyModelName::ShearViscosityConstant );
            REQUIRE( parameters.at( "modelParam" ) == 1e-5 );
         }
      }

      WHEN( "The surface tension coefficient model is read." ) {
         THEN( "The model name should be constant and the parameter should be modelParam=1e-6." ) {
            std::unordered_map<std::string, double> const parameters( reader->ReadModelData( { 1, 2 }, MaterialProperty::SurfaceTensionCoefficient ) );
            REQUIRE( reader->ReadModelName( { 1, 2 }, MaterialProperty::SurfaceTensionCoefficient ) == MaterialPropertyModelName::SurfaceTensionCoefficientConstant );
            REQUIRE( parameters.at( "modelParam" ) == 1e-6 );
         }
      }
   }

   GIVEN( "A xml document with non-existing tags to read the material data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <materials>"
                                  "  </materials>"
                                  "  <materialPairings>"
                                  "  </materialPairings>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<MaterialReader const> const reader( std::make_unique<XmlMaterialReader const>( xml_tree ) );
      WHEN( "Data is read from the tree." ) {
         THEN( "All should throw an logic error exception. The type should be fluid." ) {
            REQUIRE_THROWS_AS( reader->ReadNumberOfMaterials(), std::logic_error );
            REQUIRE( reader->ReadMaterialType( 1, MaterialType::Fluid ) == MaterialType::Fluid );
            REQUIRE_THROWS_AS( reader->ReadFixedValue( { 1 }, MaterialProperty::SpecificHeatCapacity ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadFixedValue( { 1 }, MaterialProperty::ThermalConductivity ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadFixedValue( { 1 }, MaterialProperty::BulkViscosity ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadFixedValue( { 1 }, MaterialProperty::ShearViscosity ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadFixedValue( { 1, 2 }, MaterialProperty::SurfaceTensionCoefficient ), std::logic_error );
         }
      }
   }
}
