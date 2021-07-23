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

#include <algorithm>
#include <catch2/catch.hpp>
#include <deque>
#include <iterator>
#include <list>
#include <math.h>
#include <unordered_map>
#include <vector>
#include <array>
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

SCENARIO( "Array of elementwise application is produced properly", "[1rank]" ) {
   GIVEN( "A container providing the necessary consexpr size and at functions" ) {
      constexpr std::array<double, 4> a = { 1, -2, 3, -4 };
      WHEN( "We apply two different elementwise functions" ) {
         constexpr auto square = []( auto const in ) { return in * in; };
         constexpr auto negate = []( auto const in ) { return -in; };
         //constinit would be preferred here, but lacking compiler support.
         constexpr auto array_of_a_squared = ContainerOperations::ArrayOfElementWiseFunctionApplication( a, square );
         constexpr auto array_of_a_negated = ContainerOperations::ArrayOfElementWiseFunctionApplication( a, negate );
         THEN( "The resutls match the expectation" ) {
            constexpr std::array<double, a.size()> a_squared_expected = { 1, 4, 9, 16 };
            constexpr std::array<double, a.size()> a_negated_expected = { -1, 2, -3, 4 };
            REQUIRE( array_of_a_squared == a_squared_expected );
            REQUIRE( array_of_a_negated == a_negated_expected );
         }
      }
   }
}
SCENARIO( "Transform-If transforms a container into another one if the condition is met", "[1rank]" ) {
   GIVEN( "Two containers of different type but holding same elements, and one without the common elements" ) {
      constexpr int common_element       = 42;
      std::deque<int> double_ended_queue = { 1, 2, common_element, 3, 4, common_element, 5, common_element };
      std::unordered_map<int, int> umap  = { { 1, 1 }, { 2, 2 }, { 3, common_element }, { 4, 3 }, { 5, 4 }, { 6, common_element }, { 7, 5 }, { 8, common_element } };
      std::list<int> list                = { 1, 2, 3, 4, 5 };
      WHEN( "We transform the containers into a vector under the contion to exclude the common element and add one to the remaining entries" ) {
         std::vector<int> vector_from_dequeue_transformation, vector_from_umap_transformation, vector_from_list_transformation;
         ContainerOperations::transform_if(
               std::cbegin( double_ended_queue ), std::cend( double_ended_queue ), std::back_inserter( vector_from_dequeue_transformation ), [c = common_element]( auto const& in ) { return !( in == c ); }, []( auto const& in ) { return in + 1; } );
         ContainerOperations::transform_if(
               std::cbegin( umap ), std::cend( umap ), std::back_inserter( vector_from_umap_transformation ), [c = common_element]( auto const& in ) { return !( std::get<1>( in ) == c ); }, []( auto const& in ) { return std::get<1>( in ) + 1; } );
         ContainerOperations::transform_if(
               std::begin( list ), std::end( list ), std::back_inserter( vector_from_list_transformation ), [c = common_element]( auto const& in ) { return !( in == c ); }, []( auto const& in ) { return in + 1; } );
         THEN( "The resutling vectors hold the same elements" ) {
            std::sort( std::begin( vector_from_umap_transformation ), std::end( vector_from_umap_transformation ) );
            REQUIRE( vector_from_dequeue_transformation == vector_from_umap_transformation );
            REQUIRE( vector_from_umap_transformation == vector_from_list_transformation );
         }
      }
      WHEN( "We transform the containers into a vector under the condition to exclude everything but the common element and adding one" ) {
         std::vector<int> vector_from_dequeue_transformation, vector_from_umap_transformation, vector_from_list_transformation;
         ContainerOperations::transform_if(
               std::cbegin( double_ended_queue ), std::cend( double_ended_queue ), std::back_inserter( vector_from_dequeue_transformation ), [c = common_element]( auto const& in ) { return in == c; }, []( auto const& in ) { return in + 1; } );
         ContainerOperations::transform_if(
               std::cbegin( umap ), std::cend( umap ), std::back_inserter( vector_from_umap_transformation ), [c = common_element]( auto const& in ) { return std::get<1>( in ) == c; }, []( auto const& in ) { return std::get<1>( in ) + 1; } );
         ContainerOperations::transform_if(
               std::begin( list ), std::end( list ), std::back_inserter( vector_from_list_transformation ), [c = common_element]( auto const& in ) { return in == c; }, []( auto const& in ) { return in + 1; } );
         THEN( "The vectors derived form the first two containers are identical" ) {
            std::sort( std::begin( vector_from_umap_transformation ), std::end( vector_from_umap_transformation ) );
            REQUIRE( vector_from_dequeue_transformation == vector_from_umap_transformation );
         }
         THEN( "The vector derived from the last container is empty" ) {
            REQUIRE( vector_from_list_transformation.size() == 0 );
         }
      }
   }
}
