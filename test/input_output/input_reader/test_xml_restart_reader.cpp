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

#include "input_output/input_reader/restart_reader/restart_reader.h"
#include "input_output/input_reader/restart_reader/xml_restart_reader.h"

SCENARIO( "Check that the xml restart reader works properly", "[1rank]" ) {
   GIVEN( "A xml document with the valid content to read the restart data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <restart>"
                                  "     <restore> "
                                  "       <mode> Soft </mode>"
                                  "       <fileName> file.h5 </fileName>"
                                  "     </restore>"
                                  "     <snapshots>"
                                  "       <type> Stamps </type>"
                                  "       <interval> 2631 </interval>"
                                  "       <intervalsToKeep> 4 </intervalsToKeep>"
                                  "       <stamps>"
                                  "          <ts1> 1.0 </ts1>"
                                  "          <ts2> 2.0 </ts2>"
                                  "          <ts3> 3.0 </ts3>"
                                  "       </stamps>"
                                  "     </snapshots>"
                                  "  </restart>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<RestartReader const> const reader( std::make_unique<XmlRestartReader const>( xml_tree ) );
      WHEN( "The restore mode is read from the tree." ) {
         THEN( "The restore mode should be Soft" ) {
            REQUIRE( reader->ReadRestoreMode() == RestoreMode::Soft );
         }
      }
      WHEN( "The restore filename is read from the tree." ) {
         THEN( "The restore filename should be file.h5" ) {
            REQUIRE( reader->ReadRestoreFilename() == "file.h5" );
         }
      }
      WHEN( "The snapshot type is read from the tree." ) {
         THEN( "The snapshot type should be Stamps" ) {
            REQUIRE( reader->ReadSnapshotTimesType() == SnapshotTimesType::Stamps );
         }
      }
      WHEN( "The snapshots to keep is read from the tree." ) {
         THEN( "The snapshots to keep should be 4" ) {
            REQUIRE( reader->ReadSnapshotIntervalsToKeep() == 4 );
         }
      }
      WHEN( "The snapshot interval is read from the tree." ) {
         THEN( "The snapshot interval should be 2631" ) {
            REQUIRE( reader->ReadSnapshotInterval() == 2631 );
         }
      }
      WHEN( "The snapshot stamps are read from the tree." ) {
         THEN( "The snapshot stamps should be {1,2,3}" ) {
            std::vector<double> const stamps( reader->ReadSnapshotTimeStamps() );
            REQUIRE( stamps.size() == 3 );
            REQUIRE( stamps[0] == 1.0 );
            REQUIRE( stamps[1] == 2.0 );
            REQUIRE( stamps[2] == 3.0 );
         }
      }
   }

   GIVEN( "A xml document with invalid content to read the restart data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <restart>"
                                  "     <restore> "
                                  "       <mode> Forcced </mode>"
                                  "       <fileName> </fileName>"
                                  "     </restore>"
                                  "     <snapshots>"
                                  "       <type> InvertallStamps </type>"
                                  "       <interval> -4 </interval>"
                                  "       <intervalsToKeep> -5 </intervalsToKeep>"
                                  "       <stamps>"
                                  "          <tss2> 1.0 </tss2>"
                                  "          <tss3> 3.0 </tss3>"
                                  "       </stamps>"
                                  "     </snapshots>"
                                  "  </restart>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<RestartReader const> const reader( std::make_unique<XmlRestartReader const>( xml_tree ) );
      WHEN( "Data is read from the tree." ) {
         THEN( "All should throw an exception (except time stamps, which should be empty)" ) {
            REQUIRE_THROWS_AS( reader->ReadRestoreMode(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadRestoreFilename(), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadSnapshotTimesType(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadSnapshotIntervalsToKeep(), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadSnapshotInterval(), std::invalid_argument );
            std::vector<double> const stamps( reader->ReadSnapshotTimeStamps() );
            REQUIRE( stamps.size() == 0 );
         }
      }
   }

   GIVEN( "A xml document with non-existing tags to read the restart data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <restart>"
                                  "  </restart>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<RestartReader const> const reader( std::make_unique<XmlRestartReader const>( xml_tree ) );
      WHEN( "Data is read from the tree." ) {
         THEN( "All should throw a logic error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadRestoreMode(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadRestoreFilename(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadSnapshotTimesType(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadSnapshotIntervalsToKeep(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadSnapshotInterval(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadSnapshotTimeStamps(), std::logic_error );
         }
      }
   }
}
