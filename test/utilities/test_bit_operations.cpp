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
#include "utilities/bit_operations.h"
#include <bitset>

SCENARIO( "Bit rotated works as expected", "[1rank]" ) {
   GIVEN( "Two different bitsets" ) {
      std::bitset<3> const a( 2 );
      std::bitset<7> const b( 35 );
      WHEN( "We ask for the a single rotation and for a five-fold one" ) {
         constexpr std::size_t single_rotation = 1;
         constexpr std::size_t five_rotations  = 5;
         auto const single_a                   = BitOperations::RightCircularShift( a, single_rotation );
         auto const single_b                   = BitOperations::RightCircularShift( b, single_rotation );
         auto const five_a                     = BitOperations::RightCircularShift( a, five_rotations );
         auto const five_b                     = BitOperations::RightCircularShift( b, five_rotations );
         THEN( "We obtain the expected results" ) {
            std::bitset<3> expected_a_single( 1 );
            std::bitset<3> expected_a_five( 0 );
            std::bitset<7> expected_b_single( 81 );
            std::bitset<7> expected_b_five( 13 );
            REQUIRE( single_a == expected_a_single );
            REQUIRE( five_a == expected_a_five );
            REQUIRE( single_b == expected_b_single );
            REQUIRE( five_b == expected_b_five );
         }
      }
   }
}

SCENARIO( "Sub-bitests are correctly created" ) {
   GIVEN( "Two different bitsets" ) {
      constexpr std::size_t length_a = 9;
      constexpr std::size_t length_b = 15;
      std::bitset<length_a> const a( 404 );
      std::bitset<length_b> const b( 87390 );
      WHEN( "We take the first three and last for bits as subbit" ) {
         auto a_first_three = BitOperations::SubBitset<3>( a, length_a - 3 );
         auto b_first_three = BitOperations::SubBitset<3>( b, length_b - 3 );
         auto a_last_four   = BitOperations::SubBitset<4>( a, 0 );
         auto b_last_four   = BitOperations::SubBitset<4>( b, 0 );
         THEN( "The subbits match the expectation" ) {
            std::bitset<3> expected_a_first_three( 6 );
            std::bitset<3> expected_b_first_three( 5 );
            std::bitset<4> expected_a_last_four( 4 );
            std::bitset<4> expected_b_last_four( 14 );
            REQUIRE( a_first_three == expected_a_first_three );
            REQUIRE( b_first_three == expected_b_first_three );
            REQUIRE( a_last_four == expected_a_last_four );
            REQUIRE( b_last_four == expected_b_last_four );
         }
      }
   }
}
