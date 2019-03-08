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

#include <vector>
#include <algorithm>

#include "utilities/mathematical_functions.h"

SCENARIO( "Signum function correctness", "[1rank]" ) {

   GIVEN( "Negative int" ) {
      constexpr int signed_negative_int = -256;
      WHEN( "Signum function is applied to negative signed integer" ) {
         THEN( "-1 is returned" ) {
            REQUIRE( Signum( signed_negative_int ) == -1 );
         }
      }
   }
   GIVEN( "Positive int" ) {
      constexpr int signed_positive_int = 42;
      WHEN( "Signum function is applied to positive signed integer" ) {
         THEN( "1 is returned" ) {
            REQUIRE( Signum( signed_positive_int ) == 1 );
         }
      }
   }
   GIVEN( "Zero signed int" ) {
      constexpr int signed_zero_int = 0;
      WHEN( "Signum function is applied to zero signed integer" ) {
         THEN( "0 is returned" ) {
            REQUIRE( Signum( signed_zero_int ) == 0 );
         }
      }
   }
   GIVEN( "False bool" ) {
      constexpr bool false_bool = false;
      WHEN( "Signum function is applied to false bool" ) {
         THEN( "0 is returned" ) {
            REQUIRE( Signum( false_bool ) == 0 );
         }
      }
   }
   GIVEN( "Negative double" ) {
      constexpr double signed_negative_double = -3.1;
      WHEN( "Signum function is applied to negative signed double" ) {
         THEN( "-1 is returned" ) {
            REQUIRE( Signum( signed_negative_double ) == -1 );
         }
      }
   }
   GIVEN( "Positive double" ) {
      constexpr double signed_positive_double = 1e-42;
      WHEN( "Signum function is applied to positive signed double" ) {
         THEN( "1 is returned" ) {
            REQUIRE( Signum( signed_positive_double ) == 1 );
         }
      }
   }
   GIVEN( "Zero double" ) {
      constexpr double signed_zero_double = 0;
      WHEN( "Signum function is applied to zero double" ) {
         THEN( "0 is returned" ) {
            REQUIRE( Signum( signed_zero_double ) == 0 );
         }
      }
   }
}

SCENARIO( "ConsistencyManagedSum consistency", "[1rank]" ) {

   GIVEN( "Three double values to consistently sum up" ) {
      constexpr double a = 1.0;
      constexpr double b = 1e-42;
      constexpr double c = -1.0;
      WHEN( "a, b and c are summed up using ConsistencyManagedSum" ) {
         THEN( "order of summation doesn't matter" ) {
            REQUIRE( ConsistencyManagedSum( a, b, c ) == ConsistencyManagedSum( b, c, a ) );
            REQUIRE( ConsistencyManagedSum( a, b, c ) == ConsistencyManagedSum( c, a, b ) );
            REQUIRE( ConsistencyManagedSum( a, b, c ) == ConsistencyManagedSum( b, a, c ) );
            REQUIRE( ConsistencyManagedSum( a, b, c ) == ConsistencyManagedSum( c, b, a ) );
            REQUIRE( ConsistencyManagedSum( a, b, c ) == ConsistencyManagedSum( a, c, b ) );
         }
      }
   }
}

SCENARIO( "Integer Power function (IntegerPow) consistency", "[1rank]" ) {

   GIVEN( "Base value of 2.0" ) {
      constexpr double base = 2.0;
      WHEN( "Exponent of power function is zero" ) {
         constexpr int exponent = 0;
         THEN( "Result should be 1.0" ) {
            REQUIRE( IntegerPow<exponent>( base ) == 1.0 );
         }
      }
      WHEN( "Exponent of power function is one" ) {
         constexpr int exponent = 1;
         THEN( "Result should be 2.0" ) {
            REQUIRE( IntegerPow<exponent>( base ) == 2.0 );
         }
      }
      WHEN( "Exponent of power function is 4" ) {
         constexpr int exponent = 5;
         THEN( "Result should be 32" ) {
            REQUIRE( IntegerPow<exponent>( base ) == 32.0 );
         }
      }
   }
}