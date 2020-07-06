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

#include "utilities/string_operations.h"

SCENARIO( "An indentation string with empty characters can be created", "[1rank]" ) {
   GIVEN( "A indentation width of 4" ) {
      REQUIRE( StringOperations::Indent( 4 ) == "    " );
   }
}

SCENARIO( "All spaces from strings can be removed", "[1rank]" ) {
   GIVEN( "String with single space" ) {
      std::string const single_string = "string one";
      REQUIRE( StringOperations::RemoveSpaces( single_string ) == "stringone" );
   }
   GIVEN( "String with multiple spaces" ) {
      std::string const single_string = "string first second third";
      REQUIRE( StringOperations::RemoveSpaces( single_string ) == "stringfirstsecondthird" );
   }
   GIVEN( "String with multiple contiguous spaces" ) {
      std::string const single_string = " string one  two   three    ";
      REQUIRE( StringOperations::RemoveSpaces( single_string ) == "stringonetwothree" );
   }
}

SCENARIO( "All spaces from strings can be removed and the string is converted to upper case", "[1rank]" ) {
   GIVEN( "String with single space" ) {
      std::string const single_string = "string one";
      REQUIRE( StringOperations::ToUpperCaseWithoutSpaces( single_string ) == "STRINGONE" );
   }
   GIVEN( "String with multiple spaces" ) {
      std::string const single_string = "string first second third";
      REQUIRE( StringOperations::ToUpperCaseWithoutSpaces( single_string ) == "STRINGFIRSTSECONDTHIRD" );
   }
   GIVEN( "String with multiple contiguous spaces" ) {
      std::string const single_string = " string one  two   three    ";
      REQUIRE( StringOperations::ToUpperCaseWithoutSpaces( single_string ) == "STRINGONETWOTHREE" );
   }
}

SCENARIO( "A string can be converted to a string with scientific notation", "[1rank]" ) {
   GIVEN( "Positive value cropped at two digits (floor)" ) {
      double const value = 0.12345678;
      REQUIRE( StringOperations::ToScientificNotationString( value, 2, false ) == "1.23e-01" );
   }
   GIVEN( "Positive value cropped at four digits (ceil)" ) {
      double const value = 0.12345678;
      REQUIRE( StringOperations::ToScientificNotationString( value, 2, false ) == "1.23e-01" );
   }
   GIVEN( "Positive value cropped at four digits with shifted sign space" ) {
      double const value = 0.12345678;
      REQUIRE( StringOperations::ToScientificNotationString( value, 4, true ) == " 1.2346e-01" );
   }
   GIVEN( "Negative value cropped at four digits with non shifted sign space" ) {
      double const value = -0.12345678;
      REQUIRE( StringOperations::ToScientificNotationString( value, 4, false ) == "-1.2346e-01" );
   }
   GIVEN( "Negative value cropped at four digits with shifted sign space" ) {
      double const value = -0.12345678;
      REQUIRE( StringOperations::ToScientificNotationString( value, 4, true ) == "-1.2346e-01" );
   }
}

SCENARIO( "LEading numbers from strings can be removed", "[1rank]" ) {
   GIVEN( "String with single leading number" ) {
      std::string const single_string = "1string";
      REQUIRE( StringOperations::RemoveLeadingNumbers( single_string ) == "string" );
   }
   GIVEN( "String with multiple leading numbers" ) {
      std::string const single_string = "12062string";
      REQUIRE( StringOperations::RemoveLeadingNumbers( single_string ) == "string" );
   }
   GIVEN( "String with multiple single leading number and numbers at the end" ) {
      std::string const single_string = "1string2";
      REQUIRE( StringOperations::RemoveLeadingNumbers( single_string ) == "string2" );
   }
}