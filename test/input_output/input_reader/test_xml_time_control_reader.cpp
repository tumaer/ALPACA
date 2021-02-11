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

#include "input_output/input_reader/time_control_reader/time_control_reader.h"
#include "input_output/input_reader/time_control_reader/xml_time_control_reader.h"

SCENARIO( "Check that the xml time control reader works properly", "[1rank]" ) {
   GIVEN( "A xml document with the valid content to read the time control data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <timeControl>"
                                  "     <CFLNumber> 0.6 </CFLNumber>"
                                  "     <startTime> 0.0 </startTime>"
                                  "     <endTime>   1.0 </endTime>"
                                  "  </timeControl>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<TimeControlReader const> const reader( std::make_unique<XmlTimeControlReader const>( xml_tree ) );
      WHEN( "The Cfl number is read from the tree." ) {
         THEN( "The Cfl number should be 0.6" ) {
            REQUIRE( reader->ReadCFLNumber() == 0.6 );
         }
      }
      WHEN( "The Start time is read from the tree." ) {
         THEN( "The start time should be 0.0" ) {
            REQUIRE( reader->ReadStartTime() == 0.0 );
         }
      }
      WHEN( "The End time is read from the tree." ) {
         THEN( "The end time should be 1.0" ) {
            REQUIRE( reader->ReadEndTime() == 1.0 );
         }
      }
   }
   GIVEN( "A xml document with invalid content to read the time control data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <timeControl>"
                                  "     <CFLNumber> 1.1 </CFLNumber>"
                                  "     <startTime> -1.0 </startTime>"
                                  "     <endTime> -5.0 </endTime>"
                                  "  </timeControl>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<TimeControlReader const> const reader( std::make_unique<XmlTimeControlReader const>( xml_tree ) );
      WHEN( "The Cfl number is read from the tree." ) {
         THEN( "An std::invalid_argument exception should be thrown" ) {
            REQUIRE_THROWS_AS( reader->ReadCFLNumber(), std::invalid_argument );
         }
      }
      WHEN( "The Start time is read from the tree." ) {
         THEN( "An std::invalid_argument exception should be thrown" ) {
            REQUIRE_THROWS_AS( reader->ReadStartTime(), std::invalid_argument );
         }
      }
      WHEN( "The End time is read from the tree." ) {
         THEN( "An std::invalid_argument exception should be thrown" ) {
            REQUIRE_THROWS_AS( reader->ReadEndTime(), std::invalid_argument );
         }
      }
   }

   GIVEN( "A xml document with the non-existent tags to read the time control data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <timeControl>"
                                  "  </timeControl>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<TimeControlReader const> const reader( std::make_unique<XmlTimeControlReader const>( xml_tree ) );
      WHEN( "Any data is read." ) {
         THEN( "All should throw std::logic_error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadCFLNumber(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadStartTime(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadEndTime(), std::logic_error );
         }
      }
   }
}
