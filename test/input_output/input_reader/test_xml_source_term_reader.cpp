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

#include "input_output/input_reader/source_term_reader/source_term_reader.h"
#include "input_output/input_reader/source_term_reader/xml_source_term_reader.h"

SCENARIO( "Check that the xml source term reader works properly", "[1rank]" ) {
   GIVEN( "A xml document with the valid content to read the source term data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <sourceTerms>"
                                  "     <gravity> "
                                  "       <x> 0.1 </x>"
                                  "       <y> 0.2 </y>"
                                  "       <z> -0.3 </z>"
                                  "     </gravity>"
                                  "  </sourceTerms>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<SourceTermReader const> const reader( std::make_unique<XmlSourceTermReader const>( xml_tree ) );
      WHEN( "The gravity is read from the tree." ) {
         THEN( "The gravity should be {0.1, 0.2, -0.3}" ) {
            REQUIRE( reader->ReadGravity( Direction::X ) == 0.1 );
            REQUIRE( reader->ReadGravity( Direction::Y ) == 0.2 );
            REQUIRE( reader->ReadGravity( Direction::Z ) == -0.3 );
         }
      }
   }
   GIVEN( "A xml document with invalid content to read the source term data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <sourceTerms>"
                                  "     <gravity> "
                                  "       <x> a </x>"
                                  "       <y> b </y>"
                                  "       <z> c </z>"
                                  "     </gravity>"
                                  "  </sourceTerms>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<SourceTermReader const> const reader( std::make_unique<XmlSourceTermReader const>( xml_tree ) );
      WHEN( "The gravity is read from the tree." ) {
         THEN( "All directions should thrown an std::invalid_argument" ) {
            REQUIRE_THROWS_AS( reader->ReadGravity( Direction::X ), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadGravity( Direction::Y ), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadGravity( Direction::Z ), std::invalid_argument );
         }
      }
   }

   GIVEN( "A xml document with non-existent tags to read to read the source term data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <sourceTerms>"
                                  "     <gravity> "
                                  "     </gravity>"
                                  "  </sourceTerms>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<SourceTermReader const> const reader( std::make_unique<XmlSourceTermReader const>( xml_tree ) );
      WHEN( "The gravity is read from the tree." ) {
         THEN( "All directions should thrown an std::logic_error" ) {
            REQUIRE_THROWS_AS( reader->ReadGravity( Direction::X ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadGravity( Direction::Y ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadGravity( Direction::Z ), std::logic_error );
         }
      }
   }
}
