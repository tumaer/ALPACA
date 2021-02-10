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
#ifndef ID_INFORMATION_H
#define ID_INFORMATION_H

#include <array>
#include <vector>
#include "topology/node_id_type.h"
#include "boundary_condition/boundary_specifications.h"
#include "user_specifications/compile_time_constants.h"

nid_t EastNeighborOfNodeWithId( nid_t const id );
nid_t WestNeighborOfNodeWithId( nid_t const id );
nid_t NorthNeighborOfNodeWithId( nid_t const id );
nid_t SouthNeighborOfNodeWithId( nid_t const id );
nid_t TopNeighborOfNodeWithId( nid_t const id );
nid_t BottomNeighborOfNodeWithId( nid_t const id );

nid_t GetNeighborId( nid_t const id, BoundaryLocation const location );

nid_t IdSeed();
nid_t CutHeadBit( nid_t id, unsigned int const level );
nid_t AddHeadBit( nid_t const id, unsigned int const level );

bool IsNaturalExternalBoundary( BoundaryLocation const location, nid_t const id,
                                std::array<unsigned int, 3> const number_of_nodes_on_level_zero );
bool IsExternalBoundary( BoundaryLocation const location, nid_t const id,
                         std::array<unsigned int, 3> const number_of_nodes_on_level_zero );

nid_t RawMortonIndexOfId( nid_t const id );

/**
 * @brief Gives the level that the given id is associated to. $NOT SAFE, WRONG INPUT RESULTS IN QUITE INCORRECT DATA$.
 * @param id Id of node to be evaluated.
 * @return Level of the id.
 */
constexpr unsigned int LevelOfNode( nid_t const id ) {

   nid_t working_id   = ( id >> 21 );//Cut Shadows
   unsigned int level = 0;

   while( working_id > 10 ) {
      working_id >>= 3;
      level++;
   }

   return level;
}

/**
 * @brief Gives the cell size for a node a given level.
 * @param node_size_on_level_zero The size (= size of internal cells) of a node on level zero.
 * @param level The level for which the cell size is evaluated.
 * @return The cell size of a node on the given level.
 */
constexpr double CellSizeOfLevel( double const node_size_on_level_zero, unsigned int const level ) {
   return node_size_on_level_zero / double( CC::ICX() ) / double( 1 << level );
}

/**
 * @brief Gives the cell size for the node of the given id.
 * @param id Id of the node whose cell size is to be evaluated.
 * @param node_size_on_level_zero The size (= size of internal cells) of a node on level zero.
 * @return The cell size of the node.
 */
constexpr double CellSizeOfId( nid_t const id, double const node_size_on_level_zero ) {
   return CellSizeOfLevel( node_size_on_level_zero, LevelOfNode( id ) );
}

/**
 * @brief Gives the length of the internal domain within the node of the given id.
 * @param id Id of the node whose size is to be evaluated.
 * @param node_size_on_level_zero The size (= size of internal cells) of a node on level zero.
 * @return Length of the internal domain in the node.
 */
constexpr double DomainSizeOfId( nid_t const id, double const node_size_on_level_zero ) {
   nid_t divisor = 1 << ( LevelOfNode( id ) );
   return node_size_on_level_zero / double( divisor );
}

/**
 * @brief Gives the coordinates of the coordinate system origin of a nodes internal cells, i.e. cell [CC::FICX()][CC::FICY()][CC::FICZ()].
 * @param id The id of a node for which the coordinates are calculated.
 * @param block_size The size (= size of internal cells) of a block.
 * @return The coordinates of the coordinate origin in the block, i.e. coordinates of corner of first internal cells.
 */
constexpr std::array<double, 3> DomainCoordinatesOfId( nid_t const id, double const block_size ) {
   nid_t level_operator  = id;
   unsigned int last_bit = ( level_operator & 0x1 );
   double size           = block_size;

   double x = 0;
   double y = 0;
   double z = 0;

   while( level_operator != 10 ) {
      x += last_bit * size;
      level_operator >>= 1;//Shift by one to the right
      last_bit = ( level_operator & 0x1 );
      y += last_bit * size;
      level_operator >>= 1;//Shift by one to the right
      last_bit = ( level_operator & 0x1 );
      z += last_bit * size;
      level_operator >>= 1;//Shift by one to the right
      last_bit = ( level_operator & 0x1 );
      size *= 2;
   }

   return { x, y, z };
}

/**
 * @brief Gives the unique identifer of the parent node. $NO INPUT CHECKS, CALLER MUST ENSURE CORRECT INPUT$.
 * @param id Id of the Child.
 * @return Id of the Parent.
 */
constexpr nid_t ParentIdOfNode( nid_t const id ) { return ( id >> 3 ); }

/**
 * @brief Indicates the position of the Node among its siblings. The position is encoded as integer value with 0 = bottom-south-west to 7 = top-north-east.
 * @param id The id of the node, whose position is to be determined.
 * @return The position among its siblings.
 */
constexpr unsigned int PositionOfNodeAmongSiblings( nid_t const id ) { return ( id & 0x7 ); }

/**
 * @brief Functions to indicate whether or not a node is (one of the) east (other locations respectively) most among its siblings.
 *        I.e. has (at least one) a sibling to its west (other locations respectively).
 * @param id The id ot the node.
 * @return True if it is most east (other locations respectively), False otherwise.
 */
constexpr bool EastInSiblingPack( nid_t const id ) { return ( id & 0x1 ) == 1; }
constexpr bool WestInSiblingPack( nid_t const id ) { return ( id & 0x1 ) == 0; }
constexpr bool NorthInSiblingPack( nid_t const id ) { return ( id & 0x2 ) == 2; }
constexpr bool SouthInSiblingPack( nid_t const id ) { return ( id & 0x2 ) == 0; }
constexpr bool TopInSiblingPack( nid_t const id ) { return ( id & 0x4 ) == 4; }
constexpr bool BottomInSiblingPack( nid_t const id ) { return ( id & 0x4 ) == 0; }
/**@}*/

/**
 * @brief Gives the Ids of all eight children of a Node.
 * @param id The Id of the parent node
 * @return Ids of the children in increasing order, i.e. bottom-south-west, bottom-south-east, bottom-north-west, ... ,top-north-east.
 */
inline std::vector<nid_t> IdsOfChildren( nid_t const parent_id ) {
   return { { ( parent_id << 3 ), ( parent_id << 3 ) + 1
#if DIMENSION > 1
              ,
              ( parent_id << 3 ) + 2, ( parent_id << 3 ) + 3
#endif
#if DIMENSION == 3
              ,
              ( parent_id << 3 ) + 4, ( parent_id << 3 ) + 5, ( parent_id << 3 ) + 6, ( parent_id << 3 ) + 7
#endif
   } };
}

#endif// ID_INFORMATION_H
