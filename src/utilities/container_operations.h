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

#ifndef CONTAINER_OPERATIONS_H
#define CONTAINER_OPERATIONS_H

#include <algorithm>
#include <unordered_map>
#include <utility>

namespace ContainerOperations {
   /*
    * @brief Checks if all elements in the given container are unique.
    * @param c The container holding the elements.
    * @tparam RandomAccess Container type. Must provide random access iterators.
    * @return Whether or not the elements are unique.
    * @note Copies the container argument, as it needs to be sorted and should not mess with passed-in containers.
    */
   template<typename RandomAccess>
   bool HoldsUniqueElements( RandomAccess c ) {
      std::sort( std::begin( c ), std::end( c ) );
      return std::adjacent_find( std::cbegin( c ), std::cend( c ) ) == std::cend( c );
   }

   /**
    * @brief Gives a copy of the input container rotated to the left n-times as specified.
    * @param c The original container holding the elements to be rotated and copied.
    * @param amount_of_rotations The numer of left-rotations to perform.
    * @tparam RandomAccess Container type. Must provide random access iterators.
    * @return A copy of c with rotated elements.
    * @note Undefined behavior if the number of rotations is larger than the size of the container.
    * Copies the container argument, as it needs to create a copy of all elements at some point either way.
    */
   template<typename RandomAccess>
   RandomAccess RotatedLeftCopy( RandomAccess c, std::size_t const amount_of_rotations ) {
      std::rotate( std::begin( c ), std::begin( c ) + amount_of_rotations, std::end( c ) );
      return c;
   }

   /**
    * @brief Gives the element that is present most often in the container.
    * @param non_empty The original container holding the elements.
    * @tparam Container holding the elemnts must be iterable.
    * @return The element which occurs most frequent in the container. If multiple elements occur the most, the first one (container ordering) is returned.
    * @note Function is undefined if the container is empty.
    */
   template<typename Container>
   typename Container::value_type MostFrequentElement( Container const& non_empty ) {
      std::unordered_map<typename Container::value_type, std::size_t> count;
      count.reserve( non_empty.size() );
      for( auto const& element : non_empty ) {
         count[element]++;
      }
      return std::get<0>( *std::max_element( std::cbegin( count ), std::cend( count ), []( auto const& a, auto const& b ) { return std::get<1>( a ) < std::get<1>( b ); } ) );
   }

   /**
 * @brief Applies the given function element-wise to all elements in the container and returns an array of the results.
 * @tparam Container A container holding elements providing at least a constexpr size and at function.
 * @tparam Function The function to be applied to all elements
 * @note Does not work with empty containers.
 */
   template<typename Container, typename Function>
   constexpr auto ArrayOfElementWiseFunctionApplication( Container const c, Function const f ) {
      using value_type                        = decltype( f( std::declval<typename Container::value_type>() ) );
      std::array<value_type, c.size()> result = {};
      for( std::size_t i = 0; i < c.size(); ++i ) {
         result[i] = f( c.at( i ) );
      }
      return result;
   }

}// namespace ContainerOperations

#endif//CONTAINER_OPERATIONS_H
