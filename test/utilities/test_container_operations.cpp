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
#include <deque>
#include <list>
#include <unordered_map>
#include <vector>
#include "utilities/container_operations.h"

SCENARIO( "Containers can be checked for uniqueness", "[1rank]" ) {
   GIVEN( "A unique and non-unique vector and unordered map of different types" ) {
      std::vector<double> const unique_vector  = { 1.0, 2.0, 3.0, 4.0, 5.0 };
      std::vector<int> const non_unique_vector = { 1, 2, 2, -1, 1 };
      std::deque<int> const unique_deque       = { 1, 2, 3 };
      std::deque<int> const non_unique_deque   = { 1, 2, 3, -1, 2 };
      WHEN( "We ask if the containers are unique" ) {
         bool const unique_vector_is_unique     = ContainerOperations::HoldsUniqueElements( unique_vector );
         bool const non_unique_vector_is_unique = ContainerOperations::HoldsUniqueElements( non_unique_vector );
         bool const unique_deque_is_unique      = ContainerOperations::HoldsUniqueElements( unique_deque );
         bool const non_unique_deque_is_unique  = ContainerOperations::HoldsUniqueElements( non_unique_deque );
         THEN( "We get true for the unique ones and false otherwise" ) {
            REQUIRE( unique_vector_is_unique );
            REQUIRE_FALSE( non_unique_vector_is_unique );
            REQUIRE( unique_deque_is_unique );
            REQUIRE_FALSE( non_unique_deque_is_unique );
         }
      }
   }
   GIVEN( "An empty container" ) {
      std::vector<int> const v( 0 );
      WHEN( "We ask if it is unique" ) {
         bool const is_unique = ContainerOperations::HoldsUniqueElements( v );
         THEN( "We expect it to be true" ) {
            REQUIRE( is_unique );
         }
      }
   }
}

SCENARIO( "Containers can be copy rotated left", "[1rank]" ) {
   GIVEN( "Two containers filled with some arbitrary elements" ) {
      std::vector<double> const v = { 1.0, 2.0, 3.0, 4.0 };
      std::deque<int> d           = { -1, 2, -3 };
      WHEN( "We asked for two left copy rotations" ) {
         auto const single_left_copy_rotation_v = ContainerOperations::RotatedLeftCopy( v, 1 );
         auto const three_left_copy_rotations_v = ContainerOperations::RotatedLeftCopy( v, 3 );
         auto const single_left_copy_rotation_d = ContainerOperations::RotatedLeftCopy( d, 1 );
         auto const two_left_copy_rotations_d   = ContainerOperations::RotatedLeftCopy( d, 2 );
         THEN( "The rotations match the expectation" ) {
            std::vector<double> const expected_single_rotation_v = { 2.0, 3.0, 4.0, 1.0 };
            std::vector<double> const expected_three_rotations_v = { 4.0, 1.0, 2.0, 3.0 };
            std::deque<int> const expected_single_rotation_d     = { 2, -3, -1 };
            std::deque<int> const expected_two_rotations_d       = { -3, -1, 2 };
            REQUIRE( single_left_copy_rotation_v == expected_single_rotation_v );
            REQUIRE( three_left_copy_rotations_v == expected_three_rotations_v );
            REQUIRE( single_left_copy_rotation_d == expected_single_rotation_d );
            REQUIRE( two_left_copy_rotations_d == expected_two_rotations_d );
         }
      }
   }
}

SCENARIO( "Most frequent elements can be found", "[1rank]" ) {
   constexpr int most_occuring_element = 42;
   GIVEN( "Two containers filled with some values with different occurent counts" ) {
      std::vector<int> const v        = { 1, 2, most_occuring_element, most_occuring_element, most_occuring_element, 2, 0, -3 };
      std::list<unsigned int> const l = { 1, 2, most_occuring_element, most_occuring_element, most_occuring_element };
      WHEN( "We ask for the most frequent element" ) {
         auto const most_frequent_in_v = ContainerOperations::MostFrequentElement( v );
         auto const most_frequent_in_l = ContainerOperations::MostFrequentElement( l );
         THEN( "The result is the most occuring element" ) {
            REQUIRE( most_frequent_in_v == most_occuring_element );
            REQUIRE( most_frequent_in_l == most_occuring_element );
         }
      }
   }
   GIVEN( "Two containers holding the same elements in reverse order of each other and two elemnts with highest count" ) {
      constexpr int other_most_occuring_element = 69;
      std::vector<int> const v                  = { 1, most_occuring_element, -2, most_occuring_element, 3, most_occuring_element, 4, other_most_occuring_element, -5, other_most_occuring_element, 6, other_most_occuring_element };
      std::vector<int> const w( std::crbegin( v ), std::crend( v ) );
      WHEN( "We ask for the most frequent element" ) {
         auto const most_frequent_in_v = ContainerOperations::MostFrequentElement( v );
         auto const most_frequent_in_w = ContainerOperations::MostFrequentElement( w );
         THEN( "The result is the one of the most occuring elements and the found one is opposite in the two containers" ) {
            REQUIRE( most_frequent_in_v == other_most_occuring_element );
            REQUIRE( most_frequent_in_w == most_occuring_element );
         }
      }
   }
}