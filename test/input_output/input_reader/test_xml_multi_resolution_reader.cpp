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

#include "input_output/input_reader/multi_resolution_reader/multi_resolution_reader.h"
#include "input_output/input_reader/multi_resolution_reader/xml_multi_resolution_reader.h"

SCENARIO( "Check that the xml multiresolution reader works properly", "[1rank]" ) {
   GIVEN( "A xml document with the valid content to read the restart data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <domain>"
                                  "     <nodeSize> 2.0 </nodeSize>"
                                  "     <nodeRatio>"
                                  "       <x> 5 </x>"
                                  "       <y> 7 </y>"
                                  "       <z> 2 </z>"
                                  "     </nodeRatio>"
                                  "  </domain>"
                                  "  <multiResolution>"
                                  "    <maximumLevel> 4 </maximumLevel>"
                                  "    <refinementCriterion>"
                                  "       <epsilonReference>    1e-3 </epsilonReference>"
                                  "       <levelOfEpsilonReference> 6 </levelOfEpsilonReference>"
                                  "    </refinementCriterion>"
                                  "  </multiResolution>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<MultiResolutionReader const> const reader( std::make_unique<XmlMultiResolutionReader const>( xml_tree ) );
      WHEN( "The domain data is read from the tree." ) {
         THEN( "The node size should be 2.0 and the node ratio {5,7,2}." ) {
            REQUIRE( reader->ReadNodeSizeOnLevelZero() == 2.0 );
            REQUIRE( reader->ReadNumberOfNodes( Direction::X ) == 5 );
            REQUIRE( reader->ReadNumberOfNodes( Direction::Y ) == 7 );
            REQUIRE( reader->ReadNumberOfNodes( Direction::Z ) == 2 );
         }
      }
      WHEN( "The multiresolution data is read from the tree." ) {
         THEN( "The maximum level should be 4, the epsilon reference 1e-3 and the level reference 6." ) {
            REQUIRE( reader->ReadMaximumLevel() == 4 );
            REQUIRE( reader->ReadEpsilonReference() == 1e-3 );
            REQUIRE( reader->ReadEpsilonLevelReference() == 6 );
         }
      }
   }

   GIVEN( "A xml document with invalid content to read the multiresolution data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <domain>"
                                  "     <nodeSize> -1.0 </nodeSize>"
                                  "     <nodeRatio>"
                                  "       <x> -1 </x>"
                                  "       <y> 130 </y>"
                                  "       <z> -7 </z>"
                                  "     </nodeRatio>"
                                  "  </domain>"
                                  "  <multiResolution>"
                                  "    <maximumLevel> 15 </maximumLevel>"
                                  "    <refinementCriterion>"
                                  "       <epsilonReference> -1.0 </epsilonReference>"
                                  "       <levelOfEpsilonReference> -1 </levelOfEpsilonReference>"
                                  "    </refinementCriterion>"
                                  "  </multiResolution>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<MultiResolutionReader const> const reader( std::make_unique<XmlMultiResolutionReader const>( xml_tree ) );
      WHEN( "Data is read from the tree." ) {
         THEN( "All should throw an invalid argument exception" ) {
            REQUIRE_THROWS_AS( reader->ReadNodeSizeOnLevelZero(), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadNumberOfNodes( Direction::X ), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadNumberOfNodes( Direction::Y ), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadNumberOfNodes( Direction::Z ), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadMaximumLevel(), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadEpsilonReference(), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadEpsilonLevelReference(), std::invalid_argument );
         }
      }
   }

   GIVEN( "A xml document with non-existing tags to read the multiresolution data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <domain>"
                                  "  </domain>"
                                  "  <multiResolution>"
                                  "  </multiResolution>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<MultiResolutionReader const> const reader( std::make_unique<XmlMultiResolutionReader const>( xml_tree ) );
      WHEN( "Data is read from the tree." ) {
         THEN( "All should throw an logic error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadNodeSizeOnLevelZero(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadNumberOfNodes( Direction::X ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadNumberOfNodes( Direction::Y ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadNumberOfNodes( Direction::Z ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadMaximumLevel(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadEpsilonReference(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadEpsilonLevelReference(), std::logic_error );
         }
      }
   }
}
