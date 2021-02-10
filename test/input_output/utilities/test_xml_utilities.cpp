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
#include <tinyxml2.h>

#include "input_output/utilities/xml_utilities.h"

using Catch::Matchers::Contains;

namespace {
   /**
    * @brief Creates a xml-element tree with appropriate content.
    * @return The created xml tree.
    */
   std::string CreateXmlTree() {

      return "<document>"
             "  <sec1>"
             "     <sec11>"
             "        <sec111> 1.11 </sec111>"
             "     </sec11>"
             "     <sec12> 1.2trailing chars do not matter </sec12>"
             "     <sec13> leading chars do matter1.2 </sec13>"
             "  </sec1>"
             "  <sec2>"
             "     <sec21> 9 </sec21>"
             "     <sec22> 9trail </sec22>"
             "     <sec23> lead9 </sec23>"
             "     <sec24> -2 </sec24>"
             "     <sec25> -1 </sec25>"
             "     <sec25> -2 </sec25>"
             "     <sec25> -3 </sec25>"
             "  </sec2>"
             "  <sec3>"
             "    <sec31> \n\t leading and trailing break chars \n   \t  </sec31>"
             "    <sec32>intermediate \n\t break chars</sec32>"
             "    <sec33> </sec33>"
             "  </sec3>"
             "</document>";
   }
}// namespace

SCENARIO( "Child nodes are read properly from a xml-tree", "[1rank]" ) {
   GIVEN( "A tree with seven tags on three levels." ) {
      // std::unique_ptr<tinyxml2::XMLDocument> xml_tree = CreateXmlTree();
      tinyxml2::XMLDocument xml_tree;
      xml_tree.Parse( CreateXmlTree().c_str() );
      WHEN( "The tag section 1-1-1 is queried from the tree." ) {
         THEN( "No error should be thrown." ) {
            REQUIRE_NOTHROW( XmlUtilities::GetChild( xml_tree, { "document", "sec1", "sec11", "sec111" } ) );
         }
      }
      WHEN( "The tag section 1-1-2 is queried from the tree." ) {
         THEN( "A logic_error should be thrown an the error message should contain the parent tags in angular brackets." ) {
            REQUIRE_THROWS_AS( XmlUtilities::GetChild( xml_tree, { "document", "sec1", "sec11", "sec112" } ), std::logic_error );
            REQUIRE_THROWS_WITH( XmlUtilities::GetChild( xml_tree, { "document", "sec1", "sec11", "sec112" } ),
                                 Contains( "<document>" ) && Contains( "<sec1>" ) && Contains( "<sec11>" ) && Contains( "<sec112>" ) );
         }
      }
      WHEN( "The tag section 2-1 is queried from the tree." ) {
         THEN( "No error should be thrown." ) {
            REQUIRE_NOTHROW( XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec21" } ) );
         }
      }
      WHEN( "The multiple existing tag sections 2-5 are queried from tag section 2." ) {
         tinyxml2::XMLElement const* sec2_node          = XmlUtilities::GetChild( xml_tree, { "document", "sec2" } );
         std::vector<tinyxml2::XMLElement const*> nodes = XmlUtilities::GetChilds( sec2_node, "sec25" );
         THEN( "The size of vector should be three." ) {
            REQUIRE( nodes.size() == 3 );
         }
      }
      WHEN( "The multiple non-existing tag sections 2-6 are queried from tag section 2." ) {
         tinyxml2::XMLElement const* sec2_node          = XmlUtilities::GetChild( xml_tree, { "document", "sec2" } );
         std::vector<tinyxml2::XMLElement const*> nodes = XmlUtilities::GetChilds( sec2_node, "sec26" );
         THEN( "The vector should be empty." ) {
            REQUIRE( nodes.size() == 0 );
         }
      }
      WHEN( "The tag section 1-2 is asked to exist from the document." ) {
         THEN( "The function should return true." ) {
            REQUIRE( XmlUtilities::ChildExists( xml_tree, { "document", "sec1", "sec12" } ) );
         }
      }
      WHEN( "The tag section 1-2 is asked to exist from the tag section 1." ) {
         THEN( "The function should return true." ) {
            tinyxml2::XMLElement const* sec1_node = XmlUtilities::GetChild( xml_tree, { "document", "sec1" } );
            REQUIRE( XmlUtilities::ChildExists( sec1_node, "sec12" ) );
         }
      }
      WHEN( "The tag section 2-6 is asked to exist from the document." ) {
         THEN( "The function should return false." ) {
            REQUIRE( !XmlUtilities::ChildExists( xml_tree, { "document", "sec2", "sec26" } ) );
         }
      }
      WHEN( "The tag section 2-6 is asked to exist from the tag section 2." ) {
         THEN( "The function should return false." ) {
            tinyxml2::XMLElement const* sec2_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2" } );
            REQUIRE( !XmlUtilities::ChildExists( sec2_node, "sec26" ) );
         }
      }
      WHEN( "The double value from tag section 1-1-1." ) {
         THEN( "The function should return 1.11." ) {
            tinyxml2::XMLElement const* sec111_node = XmlUtilities::GetChild( xml_tree, { "document", "sec1", "sec11", "sec111" } );
            REQUIRE( XmlUtilities::ReadDouble( sec111_node ) == 1.11 );
         }
      }
      WHEN( "The double value from tag section 1-2 with trailing characters is read." ) {
         THEN( "The function should return 1.2" ) {
            tinyxml2::XMLElement const* sec12_node = XmlUtilities::GetChild( xml_tree, { "document", "sec1", "sec12" } );
            REQUIRE( XmlUtilities::ReadDouble( sec12_node ) == 1.2 );
         }
      }
      WHEN( "The double value from tag section 1-3 with leading characters is read." ) {
         THEN( "The function should throw std::invalid_argument and the error message contains the tags." ) {
            tinyxml2::XMLElement const* sec13_node = XmlUtilities::GetChild( xml_tree, { "document", "sec1", "sec13" } );
            REQUIRE_THROWS_AS( XmlUtilities::ReadDouble( sec13_node ), std::invalid_argument );
            REQUIRE_THROWS_WITH( XmlUtilities::ReadDouble( sec13_node ),
                                 Contains( "<document>" ) && Contains( "<sec1>" ) && Contains( "<sec13>" ) && Contains( "double" ) );
         }
      }
      WHEN( "The int value from tag section 2-1 is read." ) {
         THEN( "The function should return 9." ) {
            tinyxml2::XMLElement const* sec21_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec21" } );
            REQUIRE( XmlUtilities::ReadInt( sec21_node ) == 9 );
         }
      }
      WHEN( "The int value from tag section 2-2 with trailing characters is read." ) {
         THEN( "The function should return 9." ) {
            tinyxml2::XMLElement const* sec22_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec22" } );
            REQUIRE( XmlUtilities::ReadInt( sec22_node ) == 9 );
         }
      }
      WHEN( "The int value from tag section 2-3 with leading characters is read." ) {
         THEN( "The function should throw std::invalid_argument and the error message contains the tags." ) {
            tinyxml2::XMLElement const* sec23_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec23" } );
            REQUIRE_THROWS_AS( XmlUtilities::ReadInt( sec23_node ), std::invalid_argument );
            REQUIRE_THROWS_WITH( XmlUtilities::ReadInt( sec23_node ),
                                 Contains( "<document>" ) && Contains( "<sec2>" ) && Contains( "<sec23>" ) && Contains( "int" ) );
         }
      }
      WHEN( "The int value from tag section 2-2 with trailing characters is read as long int." ) {
         THEN( "The function should return 9." ) {
            tinyxml2::XMLElement const* sec22_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec22" } );
            REQUIRE( XmlUtilities::ReadLongInt( sec22_node ) == 9 );
         }
      }
      WHEN( "The long int value from tag section 2-3 with leading characters is read." ) {
         THEN( "The function should throw std::invalid_argument and the error message contains the tags." ) {
            tinyxml2::XMLElement const* sec23_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec23" } );
            REQUIRE_THROWS_AS( XmlUtilities::ReadLongInt( sec23_node ), std::invalid_argument );
            REQUIRE_THROWS_WITH( XmlUtilities::ReadLongInt( sec23_node ),
                                 Contains( "<document>" ) && Contains( "<sec2>" ) && Contains( "<sec23>" ) && Contains( "long int" ) );
         }
      }
      WHEN( "The int value from tag section 2-1 is read as unsigned int" ) {
         THEN( "The function should return 9." ) {
            tinyxml2::XMLElement const* sec21_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec21" } );
            REQUIRE( XmlUtilities::ReadUnsignedInt( sec21_node ) == 9 );
         }
      }
      WHEN( "The negative int value from tag section 2-4 is read as unsigned int" ) {
         THEN( "The function should throw std::invalid_argument and the error message contains the tags." ) {
            tinyxml2::XMLElement const* sec24_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec24" } );
            REQUIRE_THROWS_AS( XmlUtilities::ReadUnsignedInt( sec24_node ), std::invalid_argument );
            REQUIRE_THROWS_WITH( XmlUtilities::ReadUnsignedInt( sec24_node ),
                                 Contains( "<document>" ) && Contains( "<sec2>" ) && Contains( "<sec24>" ) && Contains( "unsigned int" ) );
         }
      }
      WHEN( "The int value from tag section 2-1 is read as std::uint64_t" ) {
         THEN( "The function should return 9." ) {
            tinyxml2::XMLElement const* sec21_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec21" } );
            REQUIRE( XmlUtilities::ReadUnsignedLongInt( sec21_node ) == 9 );
         }
      }
      WHEN( "The negative int value from tag section 2-4 is read as std::uint64_t" ) {
         THEN( "The function should throw std::invalid_argument and the error message contains the tags." ) {
            tinyxml2::XMLElement const* sec24_node = XmlUtilities::GetChild( xml_tree, { "document", "sec2", "sec24" } );
            REQUIRE_THROWS_AS( XmlUtilities::ReadUnsignedLongInt( sec24_node ), std::invalid_argument );
            REQUIRE_THROWS_WITH( XmlUtilities::ReadUnsignedLongInt( sec24_node ),
                                 Contains( "<document>" ) && Contains( "<sec2>" ) && Contains( "<sec24>" ) && Contains( "unsigned long int" ) );
         }
      }
      WHEN( "The string from tag section 3-1 with leading and trailing delimiters is read." ) {
         THEN( "The string should not contain any of those delimiters." ) {
            tinyxml2::XMLElement const* sec31_node = XmlUtilities::GetChild( xml_tree, { "document", "sec3", "sec31" } );
            REQUIRE( XmlUtilities::ReadString( sec31_node ) == "leading and trailing break chars" );
         }
      }
      WHEN( "The string from tag section 3-2 with intermediate delimiters is read." ) {
         THEN( "The string should contain all of those delimiters." ) {
            tinyxml2::XMLElement const* sec32_node = XmlUtilities::GetChild( xml_tree, { "document", "sec3", "sec32" } );
            REQUIRE( XmlUtilities::ReadString( sec32_node ) == "intermediate \n\t break chars" );
         }
      }
      WHEN( "The empty string from tag section 3-3 is read." ) {
         THEN( "The function should throw std::invalid_argument and the error message contains the tags." ) {
            tinyxml2::XMLElement const* sec33_node = XmlUtilities::GetChild( xml_tree, { "document", "sec3", "sec33" } );
            REQUIRE_THROWS_AS( XmlUtilities::ReadString( sec33_node ), std::invalid_argument );
            REQUIRE_THROWS_WITH( XmlUtilities::ReadString( sec33_node ),
                                 Contains( "<document>" ) && Contains( "<sec3>" ) && Contains( "<sec33>" ) && Contains( "string" ) );
         }
      }
   }
}
