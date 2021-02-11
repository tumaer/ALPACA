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

#include "topology/space_filling_curve_index.h"
#include <algorithm>
#include <array>
#include <bitset>
#include <limits>
#include "topology/id_information.h"
#include "topology/node_id_type.h"
#include "utilities/bit_operations.h"

constexpr std::size_t bisi = 64;//bisi = Bits In Sfc Index.
static_assert( std::numeric_limits<sfcidx_t>::digits == bisi, "Space-filling curve idices type must hold 64 bits." );
static_assert( std::numeric_limits<sfcidx_t>::digits >= std::numeric_limits<nid_t>::digits, "Space-filling curve indices must not be smaller in type than node ids." );

namespace {
   /**
    * @brief Helper function to generate one long bitset by concatinating a bunch of shorter ones.
    * @param array An array holding the shorter bitsets.
    * @tparam L,M,N Sizes of the short and long bitset as well as the array.
    * @return The long bitset.
    */
   template<std::size_t L, std::size_t M, std::size_t N>
   std::bitset<L> FlattenBitsetArray( std::array<std::bitset<M>, N> const array ) {
      static_assert( L > M * N, "Cannot fit the array in a bitset of provided length" );
      std::bitset<L> flattend( 0 );
      auto const pos = []( std::size_t const i, std::size_t const j ) { return ( ( N - 1 ) * M ) - M * i + j; };
      for( std::size_t i = 0; i < array.size(); ++i ) {
         for( std::size_t j = 0; j < M; ++j ) {
            flattend[pos( i, j )] = array[i][j];
         }
      }
      return flattend;
   }
}// namespace

/**
 * @brief Generates a binary hilbert index according to \cite Butz1971 from binary coordinates, i.e. array indices.
 * @note The implemenation does not directly follow the paper but the Implementation of Doug Moore (see separate copyright) and instead of the coordinates the node ids are used.
 */
namespace BinaryHilbert {
   /* The implementation within this namespace builds upon the implementation by Doug Moore under the following copyright.
    * The present code has been strongly altered from the original version by Nils Hoppe on 29. July 2020.
    * The modifications are motivated by:
    * - reducing the library to the parts needed by Alpaca.
    * - resolving compiler warnings.
    * - Using modern std::bitset.
    * - Improving code read-ability and re-usability.
    * - Incorporations of Unit-Testability.
    */
   /*
    * hilbert.c - Computes Hilbert space-filling curve coordinates, without
    * recursion, from integer index, and vice versa, and other Hilbert-related
    * calculations.  Also known as Pi-order or Peano scan.
    *
    * Author:      Doug Moore
    *              Dept. of Computational and Applied Math
    *              Rice University
    *              http://www.caam.rice.edu/~dougm
    * Date:        Sun Feb 20 2000
    * Copyright (c) 1998-2000, Rice University
    *
    * Acknowledgement:
    * This implementation is based on the work of A. R. Butz ("Alternative
    * Algorithm for Hilbert's Space-Filling Curve", IEEE Trans. Comp., April,
    * 1971, pp 424-426) and its interpretation by Spencer W. Thomas, University
    * of Michigan (http://www-personal.umich.edu/~spencer/Home.html) in his widely
    * available C software.  While the implementation here differs considerably
    * from his, the first two interfaces and the style of some comments are very
    * much derived from his work. */
   /* LICENSE
     *
     * This software is copyrighted by Rice University.  It may be freely copied,
     * modified, and redistributed, provided that the copyright notice is
     * preserved on all copies.
     *
     * There is no warranty or other guarantee of fitness for this software,
     * it is provided solely "as is".  Bug reports or fixes may be sent
     * to the author, who may or may not act on them as he desires.
     *
     * You may include this software in a program or other software product,
     * but must display the notice:
     *
     * Hilbert Curve implementation copyright 1998, Rice University
     *
     * in any place where the end-user would see your own copyright.
     *
     * If you modify this software, you should include a notice giving the
     * name of the person performing the modification, the date of modification,
     * and the reason for such modification.
     */

   // Alpaca specific
   constexpr std::size_t group_size         = 3;
   constexpr std::size_t groups_per_node_id = 21;

   /**
    * @brief Transforms a node id into stream of bit groups as needed for the hilbert curve.
    * @param node_id The node id to be transformed.
    */
   std::bitset<bisi> GroupsStreamFromNodeId( nid_t const node_id ) {
      return std::bitset<bisi>( node_id ^ ( node_id >> group_size ) );
   }

   /**
    * @brief Gives the n-th (from left!) group of bits from a longer bitstream.
    * @param concatinated_groups The stream holding the groups contiguously.
    * @param n The group number to return.
    * @return A copy of the group.
    */
   std::bitset<group_size> NthGroupFromLeft( std::bitset<bisi> const concatinated_groups, std::size_t const n ) {
      constexpr std::size_t offset = ( groups_per_node_id - 1 ) * group_size;
      return BitOperations::SubBitset<group_size>( concatinated_groups, offset - n * group_size );
   }

   /**
    * @brief Gives the n-th bit group (from the left) from a longer bitstream after applying a flip (= xor) to it.
    * @param flip_bit_group The instruction how to flip the group.
    * @param concatinated_groups The stream holding the groups contiguously.
    * @param n The group number to return.
    * @return The flipped group.
    */
   std::bitset<group_size> NthFlippedGroup( std::bitset<group_size> const flip_bit_group, std::bitset<bisi> const concatinated_groups, std::size_t n ) {
      return flip_bit_group ^ NthGroupFromLeft( concatinated_groups, n );
   }

   /**
    * @brief Give the flip bit needed for for the given rotation.
    * @param rotation The rotation.
    * @return The corresponding flip_bit_group
    */
   std::bitset<group_size> FlipBitForRotation( std::size_t const rotation ) {
      std::bitset<group_size> flip_bit_group( 0 );
      flip_bit_group.set( rotation );
      return flip_bit_group;
   }

   /**
    * @brief Gives the rotation that should follow from the current rotation on the current group.
    * @param rotation The current rotation. Must be smaller eight.
    * @param group The current group. Group's int value mus be smaller three.
    * @return The next rotation.
    * @note Invalid inputs result in undefined behavior.
    */
   std::size_t NextRotationLookup( std::size_t const rotation, std::bitset<group_size> const group ) {
      constexpr std::array<std::array<std::size_t, 3>, 8> lookup = { { { { 1, 2, 0 } }, { { 2, 0, 1 } }, { { 0, 1, 2 } }, { { 2, 0, 1 } }, { { 1, 2, 0 } }, { { 2, 0, 1 } }, { { 0, 1, 2 } }, { { 2, 0, 1 } } } };
      return lookup[group.to_ulong()][rotation];
   }

   /**
    * @brief Gives all the groups by moving through the (adjusted) coordinates along the hilbert curve, in terms of rotating and flipping the the groups from left to right.
    * @param group_stream A stream consiting of bit groups.
    * @return an array of the correclty rotated and flipped bit groups.
    */
   std::array<std::bitset<group_size>, groups_per_node_id> DetermineGroupsFromRotations( std::bitset<bisi> const group_stream ) {
      std::array<std::bitset<group_size>, groups_per_node_id> index_groups;
      std::size_t rotation = 0;
      std::bitset<group_size> flip_bit_group( 0 );
      for( std::size_t i = 0; i < groups_per_node_id; i++ ) {
         std::bitset<group_size> const group = BitOperations::RightCircularShift( NthFlippedGroup( flip_bit_group, group_stream, i ), rotation );
         index_groups[i]                     = group;
         flip_bit_group                      = FlipBitForRotation( rotation );
         rotation                            = NextRotationLookup( rotation, group.to_ulong() );
      }
      return index_groups;
   }

   /*
    * @brief Performes an exclusive-or operation on all bit groups in the array with a constant.
    * @param groups An array holding the bit groups.
    */
   void ApplyExclusiveOrToGroups( std::array<std::bitset<group_size>, groups_per_node_id>& groups ) {
      constexpr std::bitset<group_size> xor_group( 4 );
      groups[0] ^= std::bitset<group_size>( 0 );// This line (and the start iterator shift) is obsolete as soon as node ids use 63 instead of 60 bits.
      std::for_each( std::begin( groups ) + 1, std::end( groups ), [xor_group]( std::bitset<group_size>& in ) { in ^= xor_group; } );
   }

   /*
    * @brief Core of the algorithm's logic. Uses smart xor and shifting to generate the final Hilbert Index from the rotated and flipped groups.
    * @param rotated_flipped_groups Flattened version of the adjusted bit groups.
    * @return The Hilbert Index for the given groups.
    */
   sfcidx_t IndexByShiftedSelfExclusiveOr( std::bitset<bisi> rotated_flipped_groups ) {
      for( std::size_t i = 1; i < bisi; i *= 2 )
         rotated_flipped_groups ^= std::bitset<64>( rotated_flipped_groups >> i );
      return rotated_flipped_groups.to_ulong();
   }

   /**
    * @brief Interface function to yielding the binary Hilbert index from a node id.
    * @param node_id The id of the node whose index is to be computed.
    * @return The binary Hilbert index.
    */
   sfcidx_t HilbertIndex( nid_t const node_id ) {
      std::bitset<bisi> const groups_stream                                = GroupsStreamFromNodeId( node_id );
      std::array<std::bitset<group_size>, groups_per_node_id> index_groups = DetermineGroupsFromRotations( groups_stream );
      ApplyExclusiveOrToGroups( index_groups );
      return IndexByShiftedSelfExclusiveOr( FlattenBitsetArray<bisi>( index_groups ) );
   }

}// namespace BinaryHilbert

/**
 * @brief Gives the Hilbert index for the provided node id.
 * @param node_id The id of the node whose index is to be computed.
 * @return index.
 */
sfcidx_t HilbertIndex( nid_t const node_id ) {
   return BinaryHilbert::HilbertIndex( RawMortonIndexOfId( node_id ) );
}

/**
 * @brief Gives the Lebesgue index for the provided node id.
 * @param node_id The id of the node whose index is to be computed.
 * @return index.
 */
sfcidx_t LebesgueIndex( nid_t const node_id ) {
   return RawMortonIndexOfId( node_id );
}
