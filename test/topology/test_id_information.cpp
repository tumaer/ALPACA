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
#include "topology/id_information.h"

SCENARIO( "Addition and deletion of Headbits.", "[1rank]" ) {
   GIVEN( "An Id of level = 10, inclusive of headbits" ) {
      constexpr nid_t id        = 0x50030800600000;
      constexpr nid_t morton_id = 0x30800600000;
      unsigned int const level  = LevelOfNode( id );

      WHEN( "When head bits are removed from an Id" ) {
         THEN( "Morton order of that Id is obtained" ) {
            REQUIRE( morton_id == CutHeadBit( id, level ) );
         }
      }

      WHEN( "Head bits are added to an Id" ) {
         THEN( "Headbits concatenated Morton order of that Id are obtained" ) {
            REQUIRE( id == AddHeadBit( morton_id, level ) );
         }
      }
   }
}

SCENARIO( "A 4x4x4 cubic domain's neighbor look-up at level two", "[1rank]" ) {
   GIVEN( "An Id of level = 2 whose corresponding neighbor is to be computed" ) {
      nid_t const id = AddHeadBit( 0b101010, 2 );
      WHEN( "A node on level 2 is given" ) {
         THEN( "The neighbor ids of the node in the same level are returned" ) {
            REQUIRE( GetNeighborId( id, BoundaryLocation::East ) == AddHeadBit( 0b101011, 2 ) );
            REQUIRE( GetNeighborId( id, BoundaryLocation::West ) == AddHeadBit( 0b100011, 2 ) );
            REQUIRE( GetNeighborId( id, BoundaryLocation::South ) == AddHeadBit( 0b101000, 2 ) );
            REQUIRE( GetNeighborId( id, BoundaryLocation::North ) == AddHeadBit( 0b111000, 2 ) );
            REQUIRE( GetNeighborId( id, BoundaryLocation::Top ) == AddHeadBit( 0b101110, 2 ) );
            REQUIRE( GetNeighborId( id, BoundaryLocation::Bottom ) == AddHeadBit( 0b001110, 2 ) );
         }
      }
   }
}

SCENARIO( "Neighbor lookups work", "[1rank]" ) {
   GIVEN( "A 4x4x4 cubic domain's natural external boundary and domain boundary look-up at level one" ) {
      constexpr nid_t level_two_headbit                          = 0xA000000;
      constexpr nid_t bottom_south_east_corner                   = level_two_headbit + 0b0010010;
      constexpr nid_t top_south_east_corner                      = level_two_headbit + 0b101110;
      constexpr nid_t south_west_middle_node                     = level_two_headbit + 0b000010;
      constexpr std::array<unsigned int, 3> level_zero_nodes_xyz = { 2, 2, 2 };
      WHEN( "Id of a node at external boundary is given" ) {
         THEN( "The Ids are recognized as such and not otherwise" ) {
            REQUIRE( IsNaturalExternalBoundary( BoundaryLocation::Bottom, bottom_south_east_corner, level_zero_nodes_xyz ) );
            REQUIRE( IsNaturalExternalBoundary( BoundaryLocation::Top, top_south_east_corner, level_zero_nodes_xyz ) );
         }
      }

      WHEN( "Id of a node at domain boundary is given" ) {
         THEN( "The Ids are recognized as such and not otherwise" ) {
            REQUIRE( IsExternalBoundary( BoundaryLocation::West, south_west_middle_node, level_zero_nodes_xyz ) );
            REQUIRE( IsExternalBoundary( BoundaryLocation::BottomWest, south_west_middle_node, level_zero_nodes_xyz ) );
            REQUIRE( IsExternalBoundary( BoundaryLocation::TopWest, south_west_middle_node, level_zero_nodes_xyz ) );
            REQUIRE( IsExternalBoundary( BoundaryLocation::SouthWest, south_west_middle_node, level_zero_nodes_xyz ) );
            REQUIRE( IsExternalBoundary( BoundaryLocation::NorthWest, south_west_middle_node, level_zero_nodes_xyz ) );

            REQUIRE( IsExternalBoundary( BoundaryLocation::Top, top_south_east_corner, level_zero_nodes_xyz ) );
            REQUIRE( IsExternalBoundary( BoundaryLocation::TopEast, top_south_east_corner, level_zero_nodes_xyz ) );
            REQUIRE( IsExternalBoundary( BoundaryLocation::TopWest, top_south_east_corner, level_zero_nodes_xyz ) );
         }
      }
   }
}
