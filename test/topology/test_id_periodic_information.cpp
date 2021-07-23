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
#include "topology/id_periodic_information.h"
#include "topology/id_information.h"

namespace {
   nid_t const origin_id = IdSeed();

   /**
    * @brief Gives the node ids that make up a 2x2x2 cube on level zero.
    */
   std::vector<nid_t> TwoByTwoByTwoCube() {
      return { origin_id, EastNeighborOfNodeWithId( origin_id ), NorthNeighborOfNodeWithId( origin_id ),
               EastNeighborOfNodeWithId( NorthNeighborOfNodeWithId( origin_id ) ), TopNeighborOfNodeWithId( origin_id ),
               TopNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ), TopNeighborOfNodeWithId( NorthNeighborOfNodeWithId( origin_id ) ),
               TopNeighborOfNodeWithId( NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ) ) };
   }
   /**
    * @brief Gives the node ids that make up a 4x2x1 channel on level zero.
    */
   std::vector<nid_t> FourByTwoByOneChannel() {
      return { origin_id, EastNeighborOfNodeWithId( origin_id ), EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ),
               EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ) ),
               NorthNeighborOfNodeWithId( origin_id ), NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ),
               NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ) ),
               NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( EastNeighborOfNodeWithId( origin_id ) ) ) ) };
   }
}// namespace

SCENARIO( "Neighbor lookup for 2x2x2, 4x2x1, 1x1x512 channels of different levels", "[1rank]" ) {
   GIVEN( "A 2x2x2 cube on level 0" ) {
      std::vector<nid_t> const neighbors                     = TwoByTwoByTwoCube();
      unsigned int const west_bottom_active                  = 5;
      std::array<unsigned int, 3> const level_zero_nodes_xyz = { 2, 2, 2 };
      WHEN( "The southwest corner of a 2x2x2 cube with id 000 is given" ) {
         THEN( "The neighbors of this node are east, north, northeast, top, topeast, topnorth and top northeast" ) {
            REQUIRE( neighbors[0] == AddHeadBit( 0x0, 0 ) );
            REQUIRE( neighbors[1] == AddHeadBit( 0x1, 0 ) );
            REQUIRE( neighbors[1] == GetPeriodicNeighborId( origin_id, BoundaryLocation::West, level_zero_nodes_xyz, west_bottom_active ) );
            REQUIRE( neighbors[2] == AddHeadBit( 0x2, 0 ) );
            REQUIRE( neighbors[3] == AddHeadBit( 0x3, 0 ) );
            REQUIRE( neighbors[3] == GetPeriodicNeighborId( origin_id, BoundaryLocation::NorthWest, level_zero_nodes_xyz, west_bottom_active ) );
            REQUIRE( neighbors[4] == AddHeadBit( 0x4, 0 ) );
            REQUIRE( neighbors[4] == GetPeriodicNeighborId( origin_id, BoundaryLocation::Bottom, level_zero_nodes_xyz, west_bottom_active ) );
            REQUIRE( neighbors[5] == AddHeadBit( 0x5, 0 ) );
            REQUIRE( neighbors[6] == AddHeadBit( 0x6, 0 ) );
            REQUIRE( neighbors[7] == AddHeadBit( 0x7, 0 ) );
         }
      }
   }

   GIVEN( "A 4x2x1 channel on level 0" ) {
      std::vector<nid_t> const neighbors                     = FourByTwoByOneChannel();
      unsigned int const top_west_south_active               = 7;
      std::array<unsigned int, 3> const level_zero_nodes_xyz = { 4, 2, 1 };
      WHEN( "The southwest corner of a 4x2x1 cube with id 000 is given" ) {
         THEN( "The neighbors of this node are east-east, east-east-east, north-east-east, north-east-east-east" ) {
            REQUIRE( neighbors[1] == AddHeadBit( 0x1, 0 ) );
            REQUIRE( neighbors[2] == AddHeadBit( 0x8, 0 ) );
            REQUIRE( neighbors[5] == AddHeadBit( 0x3, 0 ) );
            REQUIRE( neighbors[6] == AddHeadBit( 0xA, 0 ) );
            REQUIRE( neighbors[7] == AddHeadBit( 0xB, 0 ) );
            REQUIRE( neighbors[0] == GetPeriodicNeighborId( origin_id, BoundaryLocation::Top, level_zero_nodes_xyz, top_west_south_active ) );
            REQUIRE( neighbors[3] == GetPeriodicNeighborId( origin_id, BoundaryLocation::West, level_zero_nodes_xyz, top_west_south_active ) );
            REQUIRE( neighbors[4] == GetPeriodicNeighborId( origin_id, BoundaryLocation::South, level_zero_nodes_xyz, top_west_south_active ) );
         }
      }
   }
}
