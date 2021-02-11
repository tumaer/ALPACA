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
#include "topology/id_information.h"
#include "topology/node_id_type.h"
#include "topology/topology_manager.h"
#include "materials/material_definitions.h"
#include "communication/mpi_utilities.h"

namespace {
   nid_t const root_node_id = IdSeed();

   /**
    * @brief Refines the node on level 0 contained in a topology.
    * @param topology Instance holding the node information (indirect return).
    */
   void RefineZerothRootNode( TopologyManager& topology ) {
      topology.RefineNodeWithId( root_node_id );
      topology.UpdateTopology();
   }

   /**
    * @brief Add a single fluid to all nodes present in a topology with one single root node.
    * @param topology Instance holding the node information (indirect return).
    * @param material Material identifier that should be added.
    */
   void AddMaterialToAllNodes( TopologyManager& topology, MaterialName const material, int const rank = 0 ) {
      if( rank == 0 ) {
         auto const root_nodes = topology.GlobalIdsOnLevel( 0 );
         std::vector<nid_t> all_nodes;
         for( auto& node : root_nodes ) {
            auto&& descendants = topology.DescendantIdsOfNode( node );
            std::move( std::begin( descendants ), std::end( descendants ), std::back_inserter( all_nodes ) );
            all_nodes.push_back( node );
         }
         std::for_each( all_nodes.begin(), all_nodes.end(), [&topology, material]( nid_t const id ) { topology.AddMaterialToNode( id, material ); } );
      }
      topology.UpdateTopology();
   }

   /**
    * @brief Adds the given material to all west-most nodes on all levels in the given topology.
    * @param topology The toplogy whose nodes are to be updated.
    * @param material The material to add to the respective nodes.
    */
   void AddMaterialToWestmostNodesOnEveryLevel( TopologyManager& topology, MaterialName const material, int const rank = 0 ) {
      if( rank == 0 ) {
         std::vector<nid_t> left_most_nodes = topology.DescendantIdsOfNode( root_node_id );
         left_most_nodes.push_back( root_node_id );
         left_most_nodes.erase( std::remove_if( std::begin( left_most_nodes ),
                                                std::end( left_most_nodes ),
                                                [level_zero_nodes = topology.GetNumberOfNodesOnLevelZero()]( auto const node ) {
                                                   return !IsNaturalExternalBoundary( BoundaryLocation::West, node, level_zero_nodes );
                                                } ),
                                std::end( left_most_nodes ) );
         std::for_each( left_most_nodes.begin(), left_most_nodes.end(), [&topology, material]( nid_t const id ) { topology.AddMaterialToNode( id, material ); } );
      }
      topology.UpdateTopology();
   }
}// namespace

namespace SimplestJumpTopology {
   constexpr unsigned int node_count  = 10;
   constexpr unsigned int block_count = 10;
   constexpr unsigned int leaf_count  = 9;

   /**
    * @brief Gives (the only) parent node id.
    */
   nid_t ParentNode() {
      return root_node_id;
   }

   /**
    * @brief Gives the id of (the only) leaf pn level zero.
    */
   nid_t LevelZeroLeafNode() {
      return root_node_id + 1;
   }

   /**
    * @brief Gives the first child id.
    */
   nid_t FirstLeftChild() {
      return IdsOfChildren( root_node_id )[0];
   }

   /**
    * @brief Gives the id of the (non-existing) id of the first child of the level zero leaf.
    */
   nid_t FirstChildIdOfLeaf() {
      return IdsOfChildren( LevelZeroLeafNode() )[0];
   }

   /**
    * @brief Gives the first of the four eastmost children.
    */
   nid_t FirstEastmostChild() {
      return IdsOfChildren( root_node_id )[1];
   }

   /*
    * @brief Gives a node that has a jump on its east side.
    */
   nid_t NodeWithJumpAtEastSide() {
      return IdsOfChildren( root_node_id )[1];
   }
   /**
    * Gives all the east most children.
    */
   std::vector<nid_t> EastmostChildren() {
      return { TopNeighborOfNodeWithId( NorthNeighborOfNodeWithId( FirstEastmostChild() ) ), TopNeighborOfNodeWithId( FirstEastmostChild() ),
               NorthNeighborOfNodeWithId( FirstEastmostChild() ), FirstEastmostChild() };
   }

   namespace MaterialsAdded {
      constexpr unsigned int multinode_count     = 5;
      constexpr unsigned int interfaceleaf_count = 4;
      constexpr unsigned int block_count         = 15;
   }// namespace MaterialsAdded
   namespace Distributed {
      constexpr unsigned int nodes_blocks_on_rank_zero = 4;
      constexpr unsigned int nodes_blocks_on_rank_one  = 4;
      constexpr unsigned int nodes_blocks_on_rank_two  = 2;
      constexpr unsigned int leaves_on_rank_zero       = 4;
      constexpr unsigned int leaves_on_rank_one        = 3;
      constexpr unsigned int leaves_on_rank_two        = 2;
      constexpr int parent_rank                        = 1;
   }// namespace Distributed
   namespace MaterialsAddedDistributed {
      constexpr unsigned long long int offset_rank_zero          = 0;
      constexpr unsigned long long int node_offset_rank_one      = 5;
      constexpr unsigned long long int leaf_count                = 9;
      constexpr unsigned long long int leaf_offset_rank_one      = 5;
      constexpr unsigned long long int interface_count           = 4;
      constexpr unsigned long long int interface_offset_rank_one = 2;
      constexpr unsigned long long int block_count               = 15;
      constexpr unsigned long long int block_offset_rank_one     = 7;
   }// namespace MaterialsAddedDistributed
   namespace Parallel {
      /**
       * @brief Gives the leaf count for each of two ranks.
       * @param rank The rank to be considered.
       * @return leaf count.
       */
      constexpr std::size_t LocalLeafCount( int const rank ) {
         return rank == 0 ? 8 : rank == 1 ? 1 : 0;
      }
      /**
       * @brief Gives the node count for each of two ranks.
       * @param rank The rank to be considered.
       * @return node count.
       */
      constexpr std::size_t LocalNodeCount( int const rank ) {
         return rank == 0 ? 9 : rank == 1 ? 1 : 0;
      }

      /**
       * @brief Gives the leaf count on the given level for the given rank for each of two ranks.
       * @param rank The rank to be considered.
       * @param level The level to be considered.
       * @return leaf count.
       */
      std::size_t LocalLeafCountOnLevel( unsigned int const level, int const rank ) {
         if( rank == 0 ) {
            return level == 0 ? 0 : 8;
         } else if( rank == 1 ) {
            return level == 0 ? 1 : 0;
         } else {
            return 0;
         }
      }

      /**
       * @brief Gives the node count on the given level for the given rank for each of two ranks.
       * @param rank The rank to be considered.
       * @param level The level to be considered.
       * @return node count.
       */
      std::size_t NodeCountOnLevelOfRank( unsigned int const level, int const rank ) {
         if( rank == 0 ) {
            return level == 0 ? 1 : 8;
         } else if( rank == 1 ) {
            return level == 0 ? 1 : 0;
         } else {
            return 0;
         }
      }

      /**
       * @brief Gives the interface leaf count for the given rank for each of two ranks.
       * @param rank The rank to be considered.
       * @param level The level to be considered.
       * @return interface leaf count.
       */
      constexpr std::size_t LocalInterfaceLeafCount( int const rank ) {
         return rank == 0 ? 4 : 0;
      }
   }// namespace Parallel
}// namespace SimplestJumpTopology

SCENARIO( "Topology manager gives correct nodes/leaf/multi-phase/offset/etc. counts ", "[1rank]" ) {
   GIVEN( "A single-phase topology with eight leaves on Lmax = 1 and two nodes on level zero" ) {
      TopologyManager simplest_jump( { 2, 1, 1 }, 1 );
      RefineZerothRootNode( simplest_jump );
      AddMaterialToAllNodes( simplest_jump, MaterialName::MaterialOne );
      WHEN( "We ask for counts on this simple jump toplogy" ) {
         THEN( "We count ten nodes and nine leaves" ) {
            auto const [node_count, leaf_count] = simplest_jump.NodeAndLeafCount();
            REQUIRE( node_count == SimplestJumpTopology::node_count );
            REQUIRE( leaf_count == SimplestJumpTopology::leaf_count );
         }
         THEN( "We count ten nodes and ten blocks" ) {
            auto const [node_count, block_count] = simplest_jump.NodeAndBlockCount();
            REQUIRE( node_count == SimplestJumpTopology::node_count );
            REQUIRE( block_count == SimplestJumpTopology::block_count );
         }
      }

      WHEN( "We add a second material to some nodes" ) {
         AddMaterialToWestmostNodesOnEveryLevel( simplest_jump, MaterialName::MaterialTwo );
         THEN( "We count five multiphase nodes" ) {
            REQUIRE( simplest_jump.MultiPhaseNodeCount() == SimplestJumpTopology::MaterialsAdded::multinode_count );
         }
         THEN( "We count four interface leaves" ) {
            REQUIRE( simplest_jump.InterfaceLeafCount() == SimplestJumpTopology::MaterialsAdded::interfaceleaf_count );
         }
         THEN( "We count ten nodes and fifteen blocks" ) {
            auto const [node_count, block_count] = simplest_jump.NodeAndBlockCount();
            REQUIRE( node_count == SimplestJumpTopology::node_count );
            REQUIRE( block_count == SimplestJumpTopology::MaterialsAdded::block_count );
         }
      }

      WHEN( "We distribute the toplogy onto three ranks" ) {
         simplest_jump.PrepareLoadBalancedTopology( 3 );
         simplest_jump.UpdateTopology();
         THEN( "We count 4, 4 and 2 nodes as well as 4, 3 and 2 leaves on rank zero, one and two, respectively" ) {
            auto const nodes_and_leaves_per_rank = simplest_jump.NodesAndLeavesPerRank();
            REQUIRE( nodes_and_leaves_per_rank.size() == 3 );

            auto const [nodes_rank_zero, leaves_rank_zero] = nodes_and_leaves_per_rank[0];
            REQUIRE( nodes_rank_zero == SimplestJumpTopology::Distributed::nodes_blocks_on_rank_zero );
            REQUIRE( leaves_rank_zero == SimplestJumpTopology::Distributed::leaves_on_rank_zero );

            auto const [nodes_rank_one, leaves_rank_one] = nodes_and_leaves_per_rank[1];
            REQUIRE( nodes_rank_one == SimplestJumpTopology::Distributed::nodes_blocks_on_rank_one );
            REQUIRE( leaves_rank_one == SimplestJumpTopology::Distributed::leaves_on_rank_one );

            auto const [nodes_rank_two, leaves_rank_two] = nodes_and_leaves_per_rank[2];
            REQUIRE( nodes_rank_two == SimplestJumpTopology::Distributed::nodes_blocks_on_rank_two );
            REQUIRE( leaves_rank_two == SimplestJumpTopology::Distributed::leaves_on_rank_two );
         }
         THEN( "We count 4, 4 and 2 nodes on rank zero, one and two, respectively and the same amount of blocks" ) {
            auto const nodes_and_blocks_per_rank = simplest_jump.NodesAndBlocksPerRank();
            REQUIRE( nodes_and_blocks_per_rank.size() == 3 );

            auto const [nodes_rank_zero, blocks_rank_zero] = nodes_and_blocks_per_rank[0];
            REQUIRE( nodes_rank_zero == SimplestJumpTopology::Distributed::nodes_blocks_on_rank_zero );
            REQUIRE( blocks_rank_zero == SimplestJumpTopology::Distributed::nodes_blocks_on_rank_zero );

            auto const [nodes_rank_one, blocks_rank_one] = nodes_and_blocks_per_rank[1];
            REQUIRE( nodes_rank_one == SimplestJumpTopology::Distributed::nodes_blocks_on_rank_one );
            REQUIRE( blocks_rank_one == SimplestJumpTopology::Distributed::nodes_blocks_on_rank_one );

            auto const [nodes_rank_two, blocks_rank_two] = nodes_and_blocks_per_rank[2];
            REQUIRE( nodes_rank_two == SimplestJumpTopology::Distributed::nodes_blocks_on_rank_two );
            REQUIRE( blocks_rank_two == SimplestJumpTopology::Distributed::nodes_blocks_on_rank_two );
         }
      }

      WHEN( "We add a second material and distribute on two ranks" ) {
         AddMaterialToWestmostNodesOnEveryLevel( simplest_jump, MaterialName::MaterialTwo );
         simplest_jump.PrepareLoadBalancedTopology( 2 );
         simplest_jump.UpdateTopology();
         THEN( "We count two, and two interface leaves on rank zero and one, respectively" ) {
            auto const interface_leaves_per_rank = simplest_jump.InterfaceLeavesPerRank();
            REQUIRE( interface_leaves_per_rank.size() == 2 );
            REQUIRE( interface_leaves_per_rank[0] == 2 );
            REQUIRE( interface_leaves_per_rank[1] == 2 );
         }

         THEN( "The node offset of rank zero is zero and equal the node count on all other ranks " ) {
            REQUIRE( simplest_jump.NodeOffsetOfRank( 0 ) == 0 );
            REQUIRE( simplest_jump.NodeOffsetOfRank( 1 ) == SimplestJumpTopology::MaterialsAddedDistributed::node_offset_rank_one );
            REQUIRE( simplest_jump.NodeOffsetOfRank( 42 ) == SimplestJumpTopology::node_count );
         }

         THEN( "The leaf offset of rank zero is zero and of equal the leaf count on all others ranks" ) {
            REQUIRE( simplest_jump.LeafOffsetOfRank( 0 ) == 0 );
            REQUIRE( simplest_jump.LeafOffsetOfRank( 1 ) == SimplestJumpTopology::MaterialsAddedDistributed::leaf_offset_rank_one );
            REQUIRE( simplest_jump.LeafOffsetOfRank( 42 ) == SimplestJumpTopology::MaterialsAddedDistributed::leaf_count );
         }

         THEN( "The interface offset of rank zero is zero and equal the block count on all other ranks" ) {
            REQUIRE( simplest_jump.InterfaceLeafOffsetOfRank( 0 ) == 0 );
            REQUIRE( simplest_jump.InterfaceLeafOffsetOfRank( 1 ) == SimplestJumpTopology::MaterialsAddedDistributed::interface_offset_rank_one );
            REQUIRE( simplest_jump.InterfaceLeafOffsetOfRank( 42 ) == SimplestJumpTopology::MaterialsAddedDistributed::interface_count );
         }

         THEN( "The node/block block offset of rank zero is zero and equal the node/block count on all other ranks" ) {
            auto const [node_offset_rank_zero, block_offset_rank_zero] = simplest_jump.NodeAndBlockOffsetOfRank( 0 );
            REQUIRE( node_offset_rank_zero == SimplestJumpTopology::MaterialsAddedDistributed::offset_rank_zero );
            REQUIRE( block_offset_rank_zero == SimplestJumpTopology::MaterialsAddedDistributed::offset_rank_zero );

            auto const [node_offset_rank_one, block_offset_rank_one] = simplest_jump.NodeAndBlockOffsetOfRank( 1 );
            REQUIRE( node_offset_rank_one == SimplestJumpTopology::MaterialsAddedDistributed::node_offset_rank_one );
            REQUIRE( block_offset_rank_one == SimplestJumpTopology::MaterialsAddedDistributed::block_offset_rank_one );

            auto const [node_offset_rank_oob, block_offset_rank_oob] = simplest_jump.NodeAndBlockOffsetOfRank( 42 );
            REQUIRE( node_offset_rank_oob == SimplestJumpTopology::node_count );
            REQUIRE( block_offset_rank_oob == SimplestJumpTopology::MaterialsAddedDistributed::block_count );
         }
      }
   }
}

SCENARIO( "Single node states can be requested", "[1rank]" ) {
   GIVEN( "A single-phase topology with eight leaves on Lmax = 1 and two nodes on level zero" ) {
      TopologyManager simplest_jump( { 2, 1, 1 }, 1 );
      RefineZerothRootNode( simplest_jump );
      AddMaterialToAllNodes( simplest_jump, MaterialName::MaterialOne );
      WHEN( "We ask some node states on this toplogy" ) {
         THEN( "We correctly tell if a node exists" ) {
            nid_t const first_left_child  = SimplestJumpTopology::FirstLeftChild();
            nid_t const first_right_child = SimplestJumpTopology::FirstChildIdOfLeaf();
            REQUIRE( simplest_jump.NodeExists( first_left_child ) );
            REQUIRE_FALSE( simplest_jump.NodeExists( first_right_child ) );
         }
         THEN( "We correctly tell if a node is a leaf" ) {
            nid_t const parent = SimplestJumpTopology::ParentNode();
            nid_t const leaf   = SimplestJumpTopology::LevelZeroLeafNode();
            REQUIRE( simplest_jump.NodeIsLeaf( leaf ) );
            REQUIRE_FALSE( simplest_jump.NodeIsLeaf( parent ) );
         }
      }
      WHEN( "We distribute the topology onto three ranks" ) {
         simplest_jump.PrepareLoadBalancedTopology( 3 );
         simplest_jump.UpdateTopology();
         THEN( "We tell the parent onto the correct rank and not any other" ) {
            REQUIRE( simplest_jump.NodeIsOnRank( SimplestJumpTopology::ParentNode(), SimplestJumpTopology::Distributed::parent_rank ) );
            REQUIRE_FALSE( simplest_jump.NodeIsOnRank( SimplestJumpTopology::ParentNode(), SimplestJumpTopology::Distributed::parent_rank + 1 ) );
            REQUIRE( simplest_jump.GetRankOfNode( SimplestJumpTopology::ParentNode() ) == SimplestJumpTopology::Distributed::parent_rank );
            REQUIRE( simplest_jump.GetRankOfNode( SimplestJumpTopology::ParentNode() ) != SimplestJumpTopology::Distributed::parent_rank + 3 );
         }
      }
      WHEN( "We add another material to the parent node" ) {
         constexpr MaterialName another_material = MaterialName::MaterialTwo;
         simplest_jump.AddMaterialToNode( SimplestJumpTopology::ParentNode(), another_material );
         simplest_jump.UpdateTopology();
         THEN( "We tell the parent multi and the leaf node does not" ) {
            REQUIRE( simplest_jump.IsNodeMultiPhase( SimplestJumpTopology::ParentNode() ) );
            REQUIRE_FALSE( simplest_jump.IsNodeMultiPhase( SimplestJumpTopology::LevelZeroLeafNode() ) );
         }
         THEN( "We tell the parent contains material two and the leaf node does not" ) {
            REQUIRE( simplest_jump.NodeContainsMaterial( SimplestJumpTopology::ParentNode(), another_material ) );
            REQUIRE_FALSE( simplest_jump.NodeContainsMaterial( SimplestJumpTopology::LevelZeroLeafNode(), another_material ) );
         }
         THEN( "We find two materials in the parent node and one in the leaf" ) {
            REQUIRE( simplest_jump.GetMaterialsOfNode( SimplestJumpTopology::ParentNode() ).size() == 2 );
            REQUIRE( simplest_jump.GetMaterialsOfNode( SimplestJumpTopology::LevelZeroLeafNode() ).size() == 1 );
         }
         THEN( "Askign for a single material works only on the leaf not the parent" ) {
            REQUIRE_THROWS( simplest_jump.SingleMaterialOfNode( SimplestJumpTopology::ParentNode() ) );
            REQUIRE_NOTHROW( simplest_jump.SingleMaterialOfNode( SimplestJumpTopology::LevelZeroLeafNode() ) );
         }
      }
   }
}

SCENARIO( "Questions about the Topology can be answered", "[1rank]" ) {
   GIVEN( "A topology with eight leaves on Lmax = 1 and two nodes on level zero, which is periodic in the x-direction" ) {
      TopologyManager simplest_periodic_jump( { 2, 1, 1 }, 1, 1 );
      RefineZerothRootNode( simplest_periodic_jump );
      WHEN( "We ask some topological questions" ) {
         THEN( "We correctly tell if a face is jump" ) {
            nid_t const node_with_jump_east = SimplestJumpTopology::NodeWithJumpAtEastSide();
            REQUIRE( simplest_periodic_jump.FaceIsJump( node_with_jump_east, BoundaryLocation::East ) );
            REQUIRE_FALSE( simplest_periodic_jump.FaceIsJump( node_with_jump_east, BoundaryLocation::West ) );
         }
         THEN( "We correctly tell external boundaries" ) {
            REQUIRE_FALSE( simplest_periodic_jump.IsExternalTopologyBoundary( BoundaryLocation::West, SimplestJumpTopology::LevelZeroLeafNode() ) );
            REQUIRE( simplest_periodic_jump.IsExternalTopologyBoundary( BoundaryLocation::North, SimplestJumpTopology::LevelZeroLeafNode() ) );
         }
         THEN( "We correctly tell neighboring leaves" ) {
            REQUIRE( simplest_periodic_jump.GetNeighboringLeaves( SimplestJumpTopology::FirstLeftChild(), BoundaryLocation::East ).front() == SimplestJumpTopology::FirstEastmostChild() );
            REQUIRE( simplest_periodic_jump.GetNeighboringLeaves( SimplestJumpTopology::FirstEastmostChild(), BoundaryLocation::East ).front() == SimplestJumpTopology::LevelZeroLeafNode() );
            REQUIRE( simplest_periodic_jump.GetNeighboringLeaves( SimplestJumpTopology::LevelZeroLeafNode(), BoundaryLocation::West ) == SimplestJumpTopology::EastmostChildren() );
         }
         THEN( "We correctly tell neighboring nodes" ) {
            REQUIRE( simplest_periodic_jump.GetTopologyNeighborId( SimplestJumpTopology::ParentNode(), BoundaryLocation::East ) == SimplestJumpTopology::LevelZeroLeafNode() );
            REQUIRE( simplest_periodic_jump.GetTopologyNeighborId( SimplestJumpTopology::LevelZeroLeafNode(), BoundaryLocation::East ) == SimplestJumpTopology::ParentNode() );
         }
      }
   }
}

SCENARIO( "Overall Information about the Topology can be given", "[1rank]" ) {
   GIVEN( "A topology with eight leaves on Lmax = 2 and 3x1x2 nodes on level zero" ) {
      TopologyManager topology( { 3, 1, 2 }, 2 );
      RefineZerothRootNode( topology );
      WHEN( "We ask questions about the overall toplogy" ) {
         THEN( "The maximum level get reported as two" ) {
            REQUIRE( topology.GetMaximumLevel() == 2 );
         }
         THEN( "The number of nodes on level zero gets reported as three in x-, one in y- and two in z-direction." ) {
            auto const [x, y, z] = topology.GetNumberOfNodesOnLevelZero();
            REQUIRE( x == 3 );
            REQUIRE( y == 1 );
            REQUIRE( z == 2 );
         }
         THEN( "The the current maximum level is reported as one" ) {
            REQUIRE( topology.GetCurrentMaximumLevel() == 1 );
         }
         THEN( "Load balancing is not neccasary before" ) {
            REQUIRE_FALSE( topology.IsLoadBalancingNecessary() );
         }
      }
      WHEN( "We coarse all the children" ) {
         for( auto const node : IdsOfChildren( root_node_id ) ) {
            topology.CoarseNodeWithId( node );
         }
         THEN( "But load balancing is neccasary afterwards" ) {
            REQUIRE( topology.IsLoadBalancingNecessary() );
         }
      }
   }
}

SCENARIO( "Node listings are correctly produced", "[2rank]" ) {
   GIVEN( "A single-phase topology with eight leaves on Lmax = 1 and two nodes on level zero" ) {
      TopologyManager simplest_jump( { 2, 1, 1 }, 1 );
      RefineZerothRootNode( simplest_jump );
      int const my_rank = MpiUtilities::MyRankId();
      AddMaterialToAllNodes( simplest_jump, MaterialName::MaterialOne, my_rank );
      WHEN( "We create list on this toplogy" ) {
         constexpr unsigned int level_zero = 0;
         constexpr unsigned int level_one  = 1;
         THEN( "There are nine leafs in the leaf list" ) {
            REQUIRE( simplest_jump.LeafIds().size() == SimplestJumpTopology::leaf_count );
         }
         THEN( "Just one leaf on level zero and eight on level one are present in the leafs-per-level lists" ) {
            REQUIRE( simplest_jump.LeafIdsOnLevel( 0 ).size() == 1 );
            REQUIRE( simplest_jump.LeafIdsOnLevel( 1 ).size() == 8 );
         }
         THEN( "Eight descendants are given for the parent and none for hte level zero leaf" ) {
            REQUIRE( simplest_jump.DescendantIdsOfNode( SimplestJumpTopology::ParentNode() ).size() == 8 );
            REQUIRE( simplest_jump.DescendantIdsOfNode( SimplestJumpTopology::LevelZeroLeafNode() ).size() == 0 );
         }
         THEN( "There are two nodes on level zero and eight on level one in the global list" ) {
            REQUIRE( simplest_jump.GlobalIdsOnLevel( 0 ).size() == 2 );
            REQUIRE( simplest_jump.GlobalIdsOnLevel( 1 ).size() == 8 );
         }
         THEN( "The amount of local leaves fits on both ranks" ) {
            REQUIRE( simplest_jump.LocalLeafIds().size() == SimplestJumpTopology::Parallel::LocalLeafCount( my_rank ) );
         }
         THEN( "The amount of local leaves per level fits on both ranks" ) {
            REQUIRE( simplest_jump.LocalLeafIdsOnLevel( level_zero ).size() == SimplestJumpTopology::Parallel::LocalLeafCountOnLevel( level_zero, my_rank ) );
            REQUIRE( simplest_jump.LocalLeafIdsOnLevel( level_one ).size() == SimplestJumpTopology::Parallel::LocalLeafCountOnLevel( level_one, my_rank ) );
         }
         THEN( "The amount of local nodes fits on both ranks" ) {
            REQUIRE( simplest_jump.LocalNodeIds().size() == SimplestJumpTopology::Parallel::LocalNodeCount( my_rank ) );
         }
         THEN( "The amounts of nodes per level per rank fits on both ranks for both ranks" ) {
            constexpr int rank_zero = 0;
            constexpr int rank_one  = 1;
            REQUIRE( simplest_jump.IdsOnLevelOfRank( level_zero, rank_zero ).size() == SimplestJumpTopology::Parallel::NodeCountOnLevelOfRank( level_zero, rank_zero ) );
            REQUIRE( simplest_jump.IdsOnLevelOfRank( level_zero, rank_one ).size() == SimplestJumpTopology::Parallel::NodeCountOnLevelOfRank( level_zero, rank_one ) );
            REQUIRE( simplest_jump.IdsOnLevelOfRank( level_one, rank_zero ).size() == SimplestJumpTopology::Parallel::NodeCountOnLevelOfRank( level_one, rank_zero ) );
            REQUIRE( simplest_jump.IdsOnLevelOfRank( level_one, rank_one ).size() == SimplestJumpTopology::Parallel::NodeCountOnLevelOfRank( level_one, rank_one ) );
         }
      }
      WHEN( "We add a second material to some nodes" ) {
         AddMaterialToWestmostNodesOnEveryLevel( simplest_jump, MaterialName::MaterialTwo );
         THEN( "The amount of interface leaves fits on both ranks" ) {
            REQUIRE( simplest_jump.LocalInterfaceLeafIds().size() == SimplestJumpTopology::Parallel::LocalInterfaceLeafCount( MpiUtilities::MyRankId() ) );
         }
      }
   }
}

SCENARIO( "Topology state can be pretty formatted", "[1rank]" ) {
   GIVEN( "A topology with one root node and lmax two" ) {
      TopologyManager topology( { 2, 1, 1 }, 1 );
      WHEN( "We refine the root node once and (pretend to) distribute across three ranks" ) {
         RefineZerothRootNode( topology );
         constexpr int number_of_ranks = 3;
         topology.PrepareLoadBalancedTopology( number_of_ranks );
         topology.UpdateTopology();
         THEN( "we can get the leaf rank distribution as pretty string" ) {
            std::string const expected_pretty_string = "+++ leave rank distribution +++ Level: 0 Rank: 0 --> 1 | Rank: 1 --> 0 | Rank: 2 --> 0 |  - Level: 1 Rank: 0 --> 3 | Rank: 1 --> 3 | Rank: 2 --> 2 |  - ";
            REQUIRE( topology.LeafRankDistribution( number_of_ranks ) == expected_pretty_string );
         }
      }
   }
}

SCENARIO( "Changes in the topology are correctly processed", "[2rank]" ) {
   GIVEN( "A topology with a single root node and lmax one" ) {
      TopologyManager topology( { 1, 1, 1 }, 1 );
      int const my_rank = MpiUtilities::MyRankId();
      WHEN( "We add two material to the root node and remove the first one" ) {
         constexpr MaterialName material_one = MaterialName::MaterialOne;
         constexpr MaterialName material_two = MaterialName::MaterialTwo;
         if( my_rank == 0 ) {
            topology.AddMaterialToNode( root_node_id, material_one );
            topology.AddMaterialToNode( root_node_id, material_two );
            topology.RemoveMaterialFromNode( root_node_id, material_one );
         }
         topology.UpdateTopology();
         THEN( "We find only the second material in the node" ) {
            REQUIRE( topology.GetMaterialsOfNode( root_node_id ) == std::vector( 1, material_two ) );
         }
      }
      WHEN( "We refine the root node and then coarse the children again" ) {
         topology.RefineNodeWithId( root_node_id );
         topology.UpdateTopology();
         THEN( "We have eight leaves and nine nodes after the refinement" ) {
            REQUIRE( topology.NodeAndLeafCount() == std::pair( 9u, 8u ) );
         }
         topology.CoarseNodeWithId( root_node_id );
         THEN( "We have one leaf and one node after the coarsening" ) {
            REQUIRE( topology.NodeAndLeafCount() == std::pair( 1u, 1u ) );
         }
      }
      WHEN( "We refine the root node" ) {
         topology.RefineNodeWithId( root_node_id );
         topology.UpdateTopology();
         THEN( "All nodes are on the same rank" ) {
            REQUIRE( topology.NodesAndLeavesPerRank().size() == 1 );
            REQUIRE( topology.NodesAndLeavesPerRank().front() == std::pair( 9u, 8u ) );
         }
      }
      WHEN( "We refine the root node and load balance the topology" ) {
         topology.RefineNodeWithId( root_node_id );
         topology.UpdateTopology();
         topology.PrepareLoadBalancedTopology( 2 );
         THEN( "The leafs are evenly distributed onto the ranks" ) {
            REQUIRE( std::get<1>( topology.NodesAndLeavesPerRank().front() ) == 4 );
            REQUIRE( std::get<1>( topology.NodesAndLeavesPerRank().back() ) == 4 );
         }
      }
      WHEN( "We refine the root node, add a couple materials and create a second topology from the first one by restoring the topology" ) {
         topology.RefineNodeWithId( root_node_id );
         topology.UpdateTopology();
         topology.PrepareLoadBalancedTopology( 2 );
         AddMaterialToAllNodes( topology, MaterialName::MaterialOne, my_rank );
         AddMaterialToWestmostNodesOnEveryLevel( topology, MaterialName::MaterialTwo, my_rank );

         TopologyManager restored_topology( { 1, 1, 1 }, 1 );
         auto&& level_one_ids = topology.GlobalIdsOnLevel( 1 );
         std::vector<nid_t> ids( std::cbegin( level_one_ids ), std::cend( level_one_ids ) );
         ids.push_back( root_node_id );
         std::vector<unsigned short> number_of_materials;
         std::transform( std::begin( ids ), std::end( ids ), std::back_inserter( number_of_materials ), [&topology]( auto const id ) { return topology.GetMaterialsOfNode( id ).size(); } );
         std::vector<unsigned short> materials;
         for( auto const id : ids ) {
            for( auto const mat : topology.GetMaterialsOfNode( id ) ) {
               materials.push_back( static_cast<unsigned short>( mat ) );
            }
         }
         REQUIRE_NOTHROW( restored_topology.RestoreTopology( ids, number_of_materials, materials ) );
         THEN( "The node lists in both topologies are equal" ) {
            REQUIRE( topology.NodeAndLeafCount() == restored_topology.NodeAndLeafCount() );
         }
         THEN( "The amount of fluids in the first child is the same on both nodes" ) {
            REQUIRE( topology.GetMaterialsOfNode( IdsOfChildren( root_node_id )[0] ) == restored_topology.GetMaterialsOfNode( IdsOfChildren( root_node_id )[0] ) );
         }
      }
   }
}
