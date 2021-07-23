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
#include "topology/topology_manager.h"
#include "topology/tree.h"
#include "communication/mpi_utilities.h"
#include "topology/id_information.h"
#include "communication/internal_halo_manager.h"

SCENARIO( "Internal Halos can be updated correctly", "[1rank],[2rank]" ) {
   constexpr MaterialName material_one = MaterialName::MaterialOne;
   constexpr MaterialName material_two = MaterialName::MaterialTwo;

   GIVEN( "A single-level all multi-phase topology with 2x2x2 neighboring nodes" ) {
      constexpr unsigned int maximum_level      = 0;
      TopologyManager two_nodes_level_zero_topo = TopologyManager( { 2, 2, 2 }, maximum_level, 0 );
      Tree tree                                 = Tree( two_nodes_level_zero_topo, maximum_level, 1.0 );
      for( auto const id : two_nodes_level_zero_topo.LocalLeafIds() ) {
         two_nodes_level_zero_topo.AddMaterialToNode( id, material_one );
         two_nodes_level_zero_topo.AddMaterialToNode( id, material_two );
         tree.CreateNode( id, { material_one, material_two } );
      }
      two_nodes_level_zero_topo.UpdateTopology();

      WHEN( "Material, interface tag and levelset values are set and halo-updated" ) {
         for( auto& [id, node] : tree.FullNodeList().at( maximum_level ) ) {
            node.SetInterfaceBlock( std::make_unique<InterfaceBlock>( static_cast<double>( id ) ) );
            auto& cells          = node.GetPhaseByMaterial( material_one ).GetRightHandSideBuffer( Equation::Mass );
            auto& interface_tags = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                     cells[i][j][k]          = static_cast<double>( id );
                     interface_tags[i][j][k] = static_cast<int8_t>( PositionOfNodeAmongSiblings( id ) );
                  }
               }
            }
         }

         CommunicationManager communication = CommunicationManager( two_nodes_level_zero_topo, maximum_level );
         InternalHaloManager internal_halos = InternalHaloManager( tree, two_nodes_level_zero_topo, communication, 1 );
         internal_halos.MaterialHaloUpdateOnLevel( maximum_level, MaterialFieldType::Conservatives, false );
         internal_halos.InterfaceTagHaloUpdateOnLevel( maximum_level, InterfaceDescriptionBufferType::Reinitialized );
         internal_halos.InterfaceHaloUpdateOnLevel( maximum_level, InterfaceBlockBufferType::LevelsetRightHandSide );

         THEN( "The halo values contain the proper values from the neighbor node" ) {
            std::unordered_map<int, std::vector<BoundaryLocation>> sides_to_check_for_id( { { 0, { BoundaryLocation::East, BoundaryLocation::NorthEast, BoundaryLocation::EastNorthTop } },
                                                                                            { 1, { BoundaryLocation::West, BoundaryLocation::NorthWest, BoundaryLocation::WestNorthTop } },
                                                                                            { 2, { BoundaryLocation::East, BoundaryLocation::SouthEast, BoundaryLocation::EastSouthTop } },
                                                                                            { 3, { BoundaryLocation::West, BoundaryLocation::SouthWest, BoundaryLocation::WestSouthTop } },
                                                                                            { 4, { BoundaryLocation::East, BoundaryLocation::NorthEast, BoundaryLocation::EastNorthBottom } },
                                                                                            { 5, { BoundaryLocation::West, BoundaryLocation::NorthWest, BoundaryLocation::WestNorthBottom } },
                                                                                            { 6, { BoundaryLocation::East, BoundaryLocation::SouthEast, BoundaryLocation::EastSouthBottom } },
                                                                                            { 7, { BoundaryLocation::West, BoundaryLocation::SouthWest, BoundaryLocation::WestSouthBottom } } } );
            for( auto& [id, node] : tree.FullNodeList().at( maximum_level ) ) {
               auto const& cells          = node.GetPhaseByMaterial( material_one ).GetRightHandSideBuffer( Equation::Mass );
               auto const& interface_tags = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
               auto const& levelset       = node.GetInterfaceBlock().GetRightHandSideBuffer( InterfaceDescription::Levelset );
               for( auto const& side : sides_to_check_for_id.at( PositionOfNodeAmongSiblings( id ) ) ) {
                  auto const recv_indices = communication.GetStartIndicesHaloRecv( side );
                  auto const size         = communication.GetHaloSize( side );
                  for( int i = recv_indices[0]; i < size[0] + recv_indices[0]; ++i ) {
                     for( int j = recv_indices[1]; j < size[1] + recv_indices[1]; ++j ) {
                        for( int k = recv_indices[2]; k < size[2] + recv_indices[2]; ++k ) {
                           REQUIRE( cells[i][j][k] == static_cast<double>( GetNeighborId( id, side ) ) );
                           REQUIRE( levelset[i][j][k] == static_cast<double>( GetNeighborId( id, side ) ) );
                           REQUIRE( interface_tags[i][j][k] == static_cast<int8_t>( PositionOfNodeAmongSiblings( GetNeighborId( id, side ) ) ) );
                        }
                     }
                  }
               }
            }
         }
      }
   }

   GIVEN( "The simplest all two-phase single-jump topology" ) {
      constexpr unsigned int maximum_level = 1;
      TopologyManager simple_jump_topo     = TopologyManager( { 2, 1, 1 }, maximum_level, 0 );
      constexpr nid_t jump_parent_id       = 0x1400001;
      Tree tree                            = Tree( simple_jump_topo, maximum_level, 1.0 );
      if( simple_jump_topo.NodeIsOnRank( jump_parent_id, MpiUtilities::MyRankId() ) ) {
         simple_jump_topo.AddMaterialToNode( jump_parent_id, material_one );
         simple_jump_topo.AddMaterialToNode( jump_parent_id, material_two );
         tree.CreateNode( jump_parent_id, { material_one, material_two } );
      }
      simple_jump_topo.RefineNodeWithId( jump_parent_id );
      simple_jump_topo.UpdateTopology();
      for( auto const id : simple_jump_topo.LocalLeafIds() ) {
         simple_jump_topo.AddMaterialToNode( id, material_one );
         simple_jump_topo.AddMaterialToNode( id, material_two );
         tree.CreateNode( id, { material_one, material_two } );
      }
      simple_jump_topo.UpdateTopology();

      WHEN( "Material, interface tag and levelset values are set and halo-updated" ) {
         for( auto& [id, node] : tree.FullNodeList()[maximum_level] ) {
            node.SetInterfaceBlock( std::make_unique<InterfaceBlock>( static_cast<double>( id ) ) );
         }
         for( auto& level : tree.FullNodeList() ) {
            for( auto& [id, node] : level ) {
               auto& cells          = node.GetPhaseByMaterial( material_one ).GetRightHandSideBuffer( Equation::Mass );
               auto& interface_tags = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
               for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                  for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                     for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        cells[i][j][k]          = static_cast<double>( id );
                        interface_tags[i][j][k] = static_cast<int8_t>( PositionOfNodeAmongSiblings( id ) );
                     }
                  }
               }
            }
         }
         CommunicationManager communication = CommunicationManager( simple_jump_topo, maximum_level );
         InternalHaloManager internal_halos = InternalHaloManager( tree, simple_jump_topo, communication, 1 );
         internal_halos.MaterialHaloUpdateOnLevel( 0, MaterialFieldType::Conservatives, false );
         internal_halos.MaterialHaloUpdateOnLevel( maximum_level, MaterialFieldType::Conservatives, false );
         internal_halos.InterfaceTagHaloUpdateOnLevel( 0, InterfaceDescriptionBufferType::Reinitialized );
         internal_halos.InterfaceTagHaloUpdateOnLevel( maximum_level, InterfaceDescriptionBufferType::Reinitialized );
         internal_halos.InterfaceHaloUpdateOnLevel( maximum_level, InterfaceBlockBufferType::LevelsetRightHandSide );

         THEN( "The values in the jump halos are correct" ) {
            std::unordered_map<int, BoundaryLocation> side_to_check_for_id( {
                  { 0, BoundaryLocation::East },
                  { 1, BoundaryLocation::West },// On level 0 [and 1] (standard)
                  { 2, BoundaryLocation::West },
                  { 4, BoundaryLocation::West },
                  { 6, BoundaryLocation::West },// For jump ("other side")
                  { 3, BoundaryLocation::West },
                  { 5, BoundaryLocation::West },
                  { 7, BoundaryLocation::West }// On level 1 (standard),
            } );
            for( auto& [id, node] : tree.FullNodeList().at( maximum_level ) ) {
               auto const& cells                   = node.GetPhaseByMaterial( material_one ).GetRightHandSideBuffer( Equation::Mass );
               auto const& interface_tags          = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
               auto const& levelset                = node.GetInterfaceBlock().GetRightHandSideBuffer( InterfaceDescription::Levelset );
               auto const& side                    = side_to_check_for_id.at( PositionOfNodeAmongSiblings( id ) );
               auto const recv_indices             = communication.GetStartIndicesHaloRecv( side );
               auto const size                     = communication.GetHaloSize( side );
               bool const is_jump                  = simple_jump_topo.FaceIsJump( id, side );
               double const material_target_value  = static_cast<double>( is_jump ? GetNeighborId( ParentIdOfNode( id ), side ) : GetNeighborId( id, side ) );
               std::uint8_t const tag_target_value = is_jump ? PositionOfNodeAmongSiblings( id ) : PositionOfNodeAmongSiblings( GetNeighborId( id, side ) );
               double const levelset_target_value  = static_cast<double>( is_jump ? id : GetNeighborId( id, side ) );
               for( int i = recv_indices[0]; i < size[0] + recv_indices[0]; ++i ) {
                  for( int j = recv_indices[1]; j < size[1] + recv_indices[1]; ++j ) {
                     for( int k = recv_indices[2]; k < size[2] + recv_indices[2]; ++k ) {
                        // Approx function required since some integer IDs cannot be converted exactly to double
                        REQUIRE( cells[i][j][k] == Approx( material_target_value ) );
                        REQUIRE( interface_tags[i][j][k] == tag_target_value );
                        REQUIRE( levelset[i][j][k] == levelset_target_value );
                     }
                  }
               }
            }
         }
      }
   }
}
