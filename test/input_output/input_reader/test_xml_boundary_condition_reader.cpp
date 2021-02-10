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

#include <string>
#include <memory>

#include "boundary_condition/boundary_specifications.h"
#include "input_output/input_reader/boundary_condition_reader/boundary_condition_reader.h"
#include "input_output/input_reader/boundary_condition_reader/xml_boundary_condition_reader.h"

SCENARIO( "Check that the xml boundary condition reader works properly", "[1rank]" ) {
   GIVEN( "A xml document with the valid and invalid content to read the boundary condition data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <domain>"
                                  "     <boundaryConditions> "
                                  "       <levelSet>"
                                  "          <east> Symmetry </east> "
                                  "          <top> Wall </top> "
                                  "          <west> periodic </west> "
                                  "          <south> ZeroGradient </south> "
                                  "       </levelSet>"
                                  "       <material>"
                                  "          <east> Symmetry </east> "
                                  "          <top> Wall </top> "
                                  "          <south> ZeroGradient </south> "
                                  "          <bottom> Periodic </bottom> "
                                  "          <north> NoBoundary </north> "
                                  "          <west> FixedValue </west> "
                                  "          <valuesWest>"
                                  "             <density> 1.0 </density>"
                                  "             <pressure> 1.0 </pressure>"
                                  "             <velocityX> 1.0 </velocityX>"
                                  "             <velocityY> 1.0 </velocityY>"
                                  "             <velocityZ> 1.0 </velocityZ>"
                                  "          </valuesWest> "
                                  "       </material>"
                                  "     </boundaryConditions>"
                                  "  </domain>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<BoundaryConditionReader const> const reader( std::make_unique<XmlBoundaryConditionReader const>( xml_tree ) );
      WHEN( "The levelset boundaries are read from the tree." ) {
         THEN( "East and south should be valid, top invalid" ) {
            REQUIRE( reader->ReadLevelsetBoundaryType( BoundaryLocation::East ) == LevelSetBoundaryType::Symmetry );
            REQUIRE( reader->ReadLevelsetBoundaryType( BoundaryLocation::South ) == LevelSetBoundaryType::ZeroGradient );
            REQUIRE( reader->ReadLevelsetBoundaryType( BoundaryLocation::West ) == LevelSetBoundaryType::Periodic );
            REQUIRE_THROWS_AS( reader->ReadLevelsetBoundaryType( BoundaryLocation::Top ), std::logic_error );
         }
      }
      WHEN( "The material boundaries are read from the tree." ) {
         THEN( "All should be valid, north invalid" ) {
            REQUIRE( reader->ReadMaterialBoundaryType( BoundaryLocation::East ) == MaterialBoundaryType::Symmetry );
            REQUIRE( reader->ReadMaterialBoundaryType( BoundaryLocation::South ) == MaterialBoundaryType::ZeroGradient );
            REQUIRE( reader->ReadMaterialBoundaryType( BoundaryLocation::Bottom ) == MaterialBoundaryType::Periodic );
            REQUIRE( reader->ReadMaterialBoundaryType( BoundaryLocation::Top ) == MaterialBoundaryType::Wall );
            REQUIRE( reader->ReadMaterialBoundaryType( BoundaryLocation::West ) == MaterialBoundaryType::FixedValue );
            REQUIRE_THROWS_AS( reader->ReadMaterialBoundaryType( BoundaryLocation::North ), std::logic_error );
         }
      }
      WHEN( "The fixed material boundaries are read from the tree for the west side." ) {
         THEN( "There should be not be any exception" ) {
            REQUIRE_NOTHROW( reader->ReadMaterialFixedValueBoundaryConditions( BoundaryLocation::West ) );
         }
      }
   }

   GIVEN( "A xml document with non-existing tags to read the boundary condition data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <domain>"
                                  "     <boundaryConditions> "
                                  "     </boundaryConditions>"
                                  "  </domain>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<BoundaryConditionReader const> const reader( std::make_unique<XmlBoundaryConditionReader const>( xml_tree ) );
      WHEN( "Data is read from the tree." ) {
         THEN( "All should throw a logic error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadLevelsetBoundaryType( BoundaryLocation::North ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadMaterialBoundaryType( BoundaryLocation::South ), std::logic_error );
         }
      }
   }
}