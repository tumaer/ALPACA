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

#include "input_output/input_reader/output_reader/output_reader.h"
#include "input_output/input_reader/output_reader/xml_output_reader.h"

SCENARIO( "Check that the xml output reader works properly", "[1rank]" ) {
   GIVEN( "A xml document with the valid content to read the output data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <output>"
                                  "     <timeNamingFactor> 1e-3 </timeNamingFactor>"
                                  "     <standardOutput>"
                                  "       <type> Stamps </type>"
                                  "       <interval> 1e-4 </interval>"
                                  "       <stamps>"
                                  "          <ts1> 1.0 </ts1>"
                                  "          <ts2> 2.0 </ts2>"
                                  "          <ts3> 3.0 </ts3>"
                                  "       </stamps>"
                                  "     </standardOutput>"
                                  "     <interfaceOutput>"
                                  "       <type> StampsInterval </type>"
                                  "       <interval> 1e2 </interval>"
                                  "       <stamps>"
                                  "          <ts1> 0.1 </ts1>"
                                  "          <ts2> 2.5 </ts2>"
                                  "          <ts3> 4.5 </ts3>"
                                  "          <ts5> 5.6 </ts5>"
                                  "       </stamps>"
                                  "     </interfaceOutput>"
                                  "  </output>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<OutputReader const> const reader( std::make_unique<XmlOutputReader const>( xml_tree ) );
      WHEN( "The time naming factor is read." ) {
         THEN( "The time naming factor should be 0.001" ) {
            REQUIRE( reader->ReadTimeNamingFactor() == 0.001 );
         }
      }
      WHEN( "The standard output data is read." ) {
         THEN( "The type should be Stamps, the interval 1e-4 and the stamps {1,2,3}." ) {
            std::vector<double> const stamps( reader->ReadOutputTimeStamps( OutputType::Standard ) );
            REQUIRE( reader->ReadOutputTimesType( OutputType::Standard ) == OutputTimesType::Stamps );
            REQUIRE( reader->ReadOutputInterval( OutputType::Standard ) == 1e-4 );
            REQUIRE( stamps.size() == 3 );
            REQUIRE( stamps[0] == 1.0 );
            REQUIRE( stamps[1] == 2.0 );
            REQUIRE( stamps[2] == 3.0 );
         }
      }
      WHEN( "The debug output data is read." ) {
         THEN( "The type should be Stamps, the interval 1e-4 and the stamps {1,2,3} (same as standard)." ) {
            std::vector<double> const stamps( reader->ReadOutputTimeStamps( OutputType::Debug ) );
            REQUIRE( reader->ReadOutputTimesType( OutputType::Debug ) == OutputTimesType::Stamps );
            REQUIRE( reader->ReadOutputInterval( OutputType::Debug ) == 1e-4 );
            REQUIRE( stamps.size() == 3 );
            REQUIRE( stamps[0] == 1.0 );
            REQUIRE( stamps[1] == 2.0 );
            REQUIRE( stamps[2] == 3.0 );
         }
      }
      WHEN( "The interface output data is read." ) {
         THEN( "The type should be IntervalStamps, the interval 1e2 and the stamps {0.1,2.5,4.5}." ) {
            std::vector<double> const stamps( reader->ReadOutputTimeStamps( OutputType::Interface ) );
            REQUIRE( reader->ReadOutputTimesType( OutputType::Interface ) == OutputTimesType::IntervalStamps );
            REQUIRE( reader->ReadOutputInterval( OutputType::Interface ) == 1e2 );
            REQUIRE( stamps.size() == 3 );
            REQUIRE( stamps[0] == 0.1 );
            REQUIRE( stamps[1] == 2.5 );
            REQUIRE( stamps[2] == 4.5 );
         }
      }
   }

   GIVEN( "A xml document with invalid content to read the output data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <output>"
                                  "     <timeNamingFactor> -1e-8 </timeNamingFactor>"
                                  "     <standardOutput>"
                                  "       <type> Stammps </type>"
                                  "       <interval> -0.00001 </interval>"
                                  "       <stamps>"
                                  "          <ts1> 0.1 </ts1>"
                                  "          <ts2> 2.0 </ts2>"
                                  "          <ts3> -1.0 </ts3>"
                                  "       </stamps>"
                                  "     </standardOutput>"
                                  "     <interfaceOutput>"
                                  "       <type> Intervval </type>"
                                  "       <interval> -1e2 </interval>"
                                  "       <stamps>"
                                  "          <ts1> -2.0 </ts1>"
                                  "          <ts2> 2.5 </ts2>"
                                  "          <ts3> 4.5 </ts3>"
                                  "          <ts5> 5.6 </ts5>"
                                  "       </stamps>"
                                  "     </interfaceOutput>"
                                  "  </output>"
                                  "</configuration>" );
      // Create the xml document
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<OutputReader const> const reader( std::make_unique<XmlOutputReader const>( xml_tree ) );

      WHEN( "The standard output data is read." ) {
         THEN( "All should throw an invalid argument or logic error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadOutputTimesType( OutputType::Standard ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadOutputInterval( OutputType::Standard ), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadOutputTimeStamps( OutputType::Standard ), std::invalid_argument );
         }
      }
      WHEN( "The debug output data is read." ) {
         THEN( "All should throw an invalid argument or logic error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadOutputTimesType( OutputType::Debug ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadOutputInterval( OutputType::Debug ), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadOutputTimeStamps( OutputType::Debug ), std::invalid_argument );
         }
      }
      WHEN( "The interface output data is read." ) {
         THEN( "All should throw an invalid argument or logic error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadOutputTimesType( OutputType::Interface ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadOutputInterval( OutputType::Interface ), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadOutputTimeStamps( OutputType::Interface ), std::invalid_argument );
         }
      }
   }

   GIVEN( "A xml document with non-existing tags to read the output data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <output>"
                                  "  </output>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<OutputReader const> const reader( std::make_unique<XmlOutputReader const>( xml_tree ) );

      WHEN( "The standard output data is read." ) {
         THEN( "All should throw an logic error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadOutputTimesType( OutputType::Standard ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadOutputInterval( OutputType::Standard ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadOutputTimeStamps( OutputType::Standard ), std::logic_error );
         }
      }
      WHEN( "The debug output data is read." ) {
         THEN( "All should throw an logic error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadOutputTimesType( OutputType::Debug ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadOutputInterval( OutputType::Debug ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadOutputTimeStamps( OutputType::Debug ), std::logic_error );
         }
      }
      WHEN( "The interface output data is read." ) {
         THEN( "All should throw an logic error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadOutputTimesType( OutputType::Interface ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadOutputInterval( OutputType::Interface ), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadOutputTimeStamps( OutputType::Interface ), std::logic_error );
         }
      }
   }
}
