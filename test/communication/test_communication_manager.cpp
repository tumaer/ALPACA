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
#include "topology/topology_manager.h"
#include "communication/communication_manager.h"
#include "communication/mpi_utilities.h"

namespace ExpectedSingleJumpHaloLists {
   std::vector<std::tuple<nid_t, BoundaryLocation>> level_zero_externals = {
         { 0x1400000, BoundaryLocation::West },
         { 0x1400000, BoundaryLocation::South },
         { 0x1400000, BoundaryLocation::North },
         { 0x1400000, BoundaryLocation::Top },
         { 0x1400000, BoundaryLocation::Bottom },
         { 0x1400001, BoundaryLocation::East },
         { 0x1400001, BoundaryLocation::South },
         { 0x1400001, BoundaryLocation::North },
         { 0x1400001, BoundaryLocation::Top },
         { 0x1400001, BoundaryLocation::Bottom } };
   std::vector<std::tuple<nid_t, BoundaryLocation>> maximum_level_externals = {
         { 0xA000008, BoundaryLocation::South },
         { 0xA000008, BoundaryLocation::Bottom },
         { 0xA000009, BoundaryLocation::East },
         { 0xA000009, BoundaryLocation::South },
         { 0xA000009, BoundaryLocation::Bottom },
         { 0xA00000A, BoundaryLocation::North },
         { 0xA00000A, BoundaryLocation::Bottom },
         { 0xA00000B, BoundaryLocation::East },
         { 0xA00000B, BoundaryLocation::North },
         { 0xA00000B, BoundaryLocation::Bottom },
         { 0xA00000C, BoundaryLocation::South },
         { 0xA00000C, BoundaryLocation::Top },
         { 0xA00000D, BoundaryLocation::East },
         { 0xA00000D, BoundaryLocation::South },
         { 0xA00000D, BoundaryLocation::Top },
         { 0xA00000E, BoundaryLocation::North },
         { 0xA00000E, BoundaryLocation::Top },
         { 0xA00000F, BoundaryLocation::East },
         { 0xA00000F, BoundaryLocation::North },
         { 0xA00000F, BoundaryLocation::Top } };

   std::vector<std::tuple<nid_t, BoundaryLocation>> maximum_level_externals_rank_one_of_three = {
         { 0xA00000B, BoundaryLocation::East },
         { 0xA00000B, BoundaryLocation::North },
         { 0xA00000B, BoundaryLocation::Bottom },
         { 0xA00000E, BoundaryLocation::North },
         { 0xA00000E, BoundaryLocation::Top },
         { 0xA00000F, BoundaryLocation::East },
         { 0xA00000F, BoundaryLocation::North },
         { 0xA00000F, BoundaryLocation::Top } };

   std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> level_zero_internals = {
         { 0x1400000, BoundaryLocation::East, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0x1400001, BoundaryLocation::West, InternalBoundaryType::NoJumpBoundaryLocal },
   };
   std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> maximum_level_internals = {
         { 0xA000008, BoundaryLocation::East, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000008, BoundaryLocation::North, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000008, BoundaryLocation::Top, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000008, BoundaryLocation::NorthEast, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000008, BoundaryLocation::TopEast, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000008, BoundaryLocation::TopNorth, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000008, BoundaryLocation::EastNorthTop, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000009, BoundaryLocation::West, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000009, BoundaryLocation::North, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000009, BoundaryLocation::Top, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000009, BoundaryLocation::NorthWest, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000009, BoundaryLocation::TopWest, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000009, BoundaryLocation::TopNorth, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA000009, BoundaryLocation::WestNorthTop, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::East, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::Top, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::SouthEast, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::TopEast, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::TopSouth, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::EastSouthTop, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000B, BoundaryLocation::West, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000B, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000B, BoundaryLocation::Top, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000B, BoundaryLocation::SouthWest, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000B, BoundaryLocation::TopWest, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000B, BoundaryLocation::TopSouth, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000B, BoundaryLocation::WestSouthTop, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::East, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::North, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::Bottom, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::NorthEast, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::BottomEast, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::BottomNorth, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::EastNorthBottom, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000D, BoundaryLocation::West, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000D, BoundaryLocation::North, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000D, BoundaryLocation::Bottom, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000D, BoundaryLocation::NorthWest, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000D, BoundaryLocation::BottomWest, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000D, BoundaryLocation::BottomNorth, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000D, BoundaryLocation::WestNorthBottom, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::East, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::Bottom, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::SouthEast, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::BottomEast, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::BottomSouth, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::EastSouthBottom, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000F, BoundaryLocation::West, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000F, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000F, BoundaryLocation::Bottom, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000F, BoundaryLocation::SouthWest, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000F, BoundaryLocation::BottomWest, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000F, BoundaryLocation::BottomSouth, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000F, BoundaryLocation::WestSouthBottom, InternalBoundaryType::NoJumpBoundaryLocal } };

   std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> maximum_level_internals_rank_one_of_three = {
         { 0xA00000B, BoundaryLocation::Top, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000B, BoundaryLocation::TopWest, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::East, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::BottomEast, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000F, BoundaryLocation::West, InternalBoundaryType::NoJumpBoundaryLocal },
         { 0xA00000F, BoundaryLocation::Bottom, InternalBoundaryType::NoJumpBoundaryLocal } };

   std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> maximum_level_internals_mpi_rank_one_of_three = {
         { 0xA00000B, BoundaryLocation::West, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000B, BoundaryLocation::West, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000B, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000B, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000B, BoundaryLocation::TopSouth, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000B, BoundaryLocation::TopSouth, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000B, BoundaryLocation::SouthWest, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000B, BoundaryLocation::SouthWest, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000B, BoundaryLocation::WestSouthTop, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000B, BoundaryLocation::WestSouthTop, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000E, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000E, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000E, BoundaryLocation::Bottom, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000E, BoundaryLocation::Bottom, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000E, BoundaryLocation::BottomSouth, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000E, BoundaryLocation::BottomSouth, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000E, BoundaryLocation::SouthEast, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000E, BoundaryLocation::SouthEast, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000E, BoundaryLocation::EastSouthBottom, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000E, BoundaryLocation::EastSouthBottom, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000F, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000F, BoundaryLocation::South, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000F, BoundaryLocation::BottomSouth, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000F, BoundaryLocation::BottomSouth, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000F, BoundaryLocation::BottomWest, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000F, BoundaryLocation::BottomWest, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000F, BoundaryLocation::SouthWest, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000F, BoundaryLocation::SouthWest, InternalBoundaryType::NoJumpBoundaryMpiRecv },
         { 0xA00000F, BoundaryLocation::WestSouthBottom, InternalBoundaryType::NoJumpBoundaryMpiSend },
         { 0xA00000F, BoundaryLocation::WestSouthBottom, InternalBoundaryType::NoJumpBoundaryMpiRecv } };

   std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> maximum_level_jumps = {
         { 0xA000008, BoundaryLocation::West, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA000008, BoundaryLocation::NorthWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA000008, BoundaryLocation::TopWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA000008, BoundaryLocation::WestNorthTop, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::West, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::SouthWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::TopWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000A, BoundaryLocation::WestSouthTop, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::West, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::NorthWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::BottomWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000C, BoundaryLocation::WestNorthBottom, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::West, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::SouthWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::BottomWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::WestSouthBottom, InternalBoundaryType::JumpBoundaryLocal } };

   std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> maximum_level_jumps_mpi_rank_one_of_three = {
         { 0xA000008, BoundaryLocation::West, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA000008, BoundaryLocation::TopWest, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA000008, BoundaryLocation::NorthWest, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA000008, BoundaryLocation::WestNorthTop, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA00000A, BoundaryLocation::West, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA00000A, BoundaryLocation::TopWest, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA00000A, BoundaryLocation::SouthWest, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA00000A, BoundaryLocation::WestSouthTop, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA00000C, BoundaryLocation::West, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA00000C, BoundaryLocation::BottomWest, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA00000C, BoundaryLocation::NorthWest, InternalBoundaryType::JumpBoundaryMpiSend },
         { 0xA00000C, BoundaryLocation::WestNorthBottom, InternalBoundaryType::JumpBoundaryMpiSend } };

   std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> maximum_level_jumps_rank_one_of_three = {
         { 0xA00000E, BoundaryLocation::West, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::BottomWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::SouthWest, InternalBoundaryType::JumpBoundaryLocal },
         { 0xA00000E, BoundaryLocation::WestSouthBottom, InternalBoundaryType::JumpBoundaryLocal } };
}// namespace ExpectedSingleJumpHaloLists

namespace {
   template<typename T>
   void RequireVectorEquality( std::vector<T> a, std::vector<T>& b ) {
      REQUIRE( a.size() == b.size() );
      std::sort( std::begin( a ), std::end( a ) );
      std::sort( std::begin( b ), std::end( b ) );
      REQUIRE( a == b );
   }
}// namespace

SCENARIO( "Communication chache is properly (in-)validated", "[1rank],[2rank]" ) {
   GIVEN( "The simplest two level communication manager" ) {
      constexpr unsigned int maximum_level = 1;
      TopologyManager topology             = TopologyManager();
      CommunicationManager communication   = CommunicationManager( topology, maximum_level );
      WHEN( "No further action is taken" ) {
         THEN( "The boundary-relation cache is invalid on all level" ) {
            REQUIRE_FALSE( communication.AreBoundariesValid( 0 ) );
            REQUIRE_FALSE( communication.AreBoundariesValid( maximum_level ) );
         }
      }
      WHEN( "The boundary relations are computed on level zero and then the cache is invalidated" ) {
         communication.GenerateNeighborRelationForHaloUpdate( 0 );
         THEN( "The boundary-relation cache is valid only on level zero" ) {
            REQUIRE( communication.AreBoundariesValid( 0 ) );
            REQUIRE_FALSE( communication.AreBoundariesValid( maximum_level ) );
         }
         communication.InvalidateCache();
         THEN( "The boundary-relation cache is invalid on all level" ) {
            REQUIRE_FALSE( communication.AreBoundariesValid( 0 ) );
            REQUIRE_FALSE( communication.AreBoundariesValid( maximum_level ) );
         }
      }
      WHEN( "The boundary relations are computed on all levels then the cache is invalidated" ) {
         communication.GenerateNeighborRelationForHaloUpdate( 0 );
         communication.GenerateNeighborRelationForHaloUpdate( maximum_level );
         THEN( "The boundary-relation cache is valid on all levels" ) {
            REQUIRE( communication.AreBoundariesValid( 0 ) );
            REQUIRE( communication.AreBoundariesValid( maximum_level ) );
         }
         communication.InvalidateCache();
         THEN( "The boundary-relation cache is invalid on all levels" ) {
            REQUIRE_FALSE( communication.AreBoundariesValid( 0 ) );
            REQUIRE_FALSE( communication.AreBoundariesValid( maximum_level ) );
         }
      }
   }
}

SCENARIO( "The neighborhood relations are correct", "[1rank]" ) {

   constexpr MaterialName material = MaterialName::MaterialOne;

   GIVEN( "A topology with two nodes in x-direction, one refined node (single jump boundary) and an empty tree" ) {
      constexpr unsigned int level_zero    = 0;
      constexpr unsigned int maximum_level = 1;
      TopologyManager simplest_jump_topo   = TopologyManager( { 2, 1, 1 }, maximum_level, 0 );
      simplest_jump_topo.AddMaterialToNode( 0x1400001, material );
      simplest_jump_topo.RefineNodeWithId( 0x1400001 );
      simplest_jump_topo.UpdateTopology();
      for( auto const& id : simplest_jump_topo.LocalLeafIds() ) {
         simplest_jump_topo.AddMaterialToNode( id, material );
      }
      simplest_jump_topo.UpdateTopology();
      WHEN( "Setting up the communication manager and generating the neighbor relations" ) {
         CommunicationManager communication = CommunicationManager( simplest_jump_topo, maximum_level );
         communication.GenerateNeighborRelationForHaloUpdate( level_zero );
         communication.GenerateNeighborRelationForHaloUpdate( maximum_level );
         THEN( "The external boundaries on are correct" ) {
            RequireVectorEquality( communication.ExternalBoundaries( level_zero ), ExpectedSingleJumpHaloLists::level_zero_externals );
            RequireVectorEquality( communication.ExternalBoundaries( maximum_level ), ExpectedSingleJumpHaloLists::maximum_level_externals );
         }
         THEN( "The internal no-jump boundaries are correct" ) {
            RequireVectorEquality( communication.InternalBoundaries( level_zero ), ExpectedSingleJumpHaloLists::level_zero_internals );
            RequireVectorEquality( communication.InternalBoundaries( maximum_level ), ExpectedSingleJumpHaloLists::maximum_level_internals );
         }
         THEN( "The internal jump boundaries are correct" ) {
            RequireVectorEquality( communication.InternalBoundariesJump( maximum_level ), ExpectedSingleJumpHaloLists::maximum_level_jumps );
         }
      }
      WHEN( "We pretend to have a three-rank topology" ) {
         simplest_jump_topo.GetLoadBalancedTopology( 3 );
         simplest_jump_topo.UpdateTopology();
         CommunicationManager communication = CommunicationManager( simplest_jump_topo, maximum_level );
         communication.GenerateNeighborRelationForHaloUpdate( level_zero );
         communication.GenerateNeighborRelationForHaloUpdate( maximum_level );
         THEN( "The external boundaries are correct" ) {
            RequireVectorEquality( communication.ExternalBoundaries( level_zero ), ExpectedSingleJumpHaloLists::level_zero_externals );
            RequireVectorEquality( communication.ExternalBoundaries( maximum_level ), ExpectedSingleJumpHaloLists::maximum_level_externals_rank_one_of_three );
         }
         THEN( "The internal no-jumps are correctly sorted into MPI and non-MPI" ) {
            REQUIRE( communication.InternalBoundariesMpi( level_zero ).size() == 0 );
            RequireVectorEquality( communication.InternalBoundaries( level_zero ), ExpectedSingleJumpHaloLists::level_zero_internals );

            RequireVectorEquality( communication.InternalBoundariesMpi( maximum_level ), ExpectedSingleJumpHaloLists::maximum_level_internals_mpi_rank_one_of_three );
            RequireVectorEquality( communication.InternalBoundaries( maximum_level ), ExpectedSingleJumpHaloLists::maximum_level_internals_rank_one_of_three );
         }
         THEN( "The internal jumps are correctly sorted into MPI and non-MPI" ) {
            RequireVectorEquality( communication.InternalBoundariesJump( maximum_level ), ExpectedSingleJumpHaloLists::maximum_level_jumps_rank_one_of_three );
            RequireVectorEquality( communication.InternalBoundariesJumpMpi( maximum_level ), ExpectedSingleJumpHaloLists::maximum_level_jumps_mpi_rank_one_of_three );
         }
      }
   }
}

SCENARIO( "Testing send and receive functions in multi-core settings", "[2rank]" ) {

   GIVEN( "A(ny) topology, a(ny) tree and a communication manager" ) {
      constexpr unsigned int maximum_level = 0;
      TopologyManager topology             = TopologyManager( { 1, 1, 1 }, maximum_level, 0 );
      CommunicationManager communicator    = CommunicationManager( topology, maximum_level );
      WHEN( "We send and receive some arbitrary data according to the communicator function descriptions" ) {
         int const my_rank      = MpiUtilities::MyRankId();
         constexpr double value = 42.0;
         std::vector<double> const data_to_send( 2, value );
         std::vector<double> received_data( 2, -1.0 );
         REQUIRE( data_to_send.size() == received_data.size() );
         std::vector<MPI_Request> requests;
         // communicators send/recv function need to be called in the correct order, otherwise tags don't match
         if( my_rank == 0 ) {
            communicator.Send( data_to_send.data(), data_to_send.size(), MPI_DOUBLE, 1, requests );
            communicator.Recv( received_data.data(), received_data.size(), MPI_DOUBLE, 1, requests );
         } else {
            communicator.Recv( received_data.data(), received_data.size(), MPI_DOUBLE, 0, requests );
            communicator.Send( data_to_send.data(), data_to_send.size(), MPI_DOUBLE, 0, requests );
         }
         MPI_Waitall( requests.size(), requests.data(), MPI_STATUSES_IGNORE );
         THEN( "The receiving vector holds the sent values" ) {
            REQUIRE( received_data.at( 0 ) == value );
            REQUIRE( received_data.at( 1 ) == value );
         }
      }
   }
}