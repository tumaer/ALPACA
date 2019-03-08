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
#include "materials/material_definitions.h"

namespace {
   constexpr std::uint64_t root_node_id = 0x1400000;

   /**
    * @brief Generates a topology with one node on level 0.
    * @param maximum_level Maximum number of levels desired.
    * @return Created TopologyManager.
    */
   TopologyManager SingleRootNodeTopoloyWithMaximumLevel( unsigned int const maximum_level = 1 ) {
      return TopologyManager( { 1, 1, 1 }, maximum_level, 0 );
   }
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
   void AddMaterialToAllNodes( TopologyManager& topology, MaterialName const material ) {
      auto all_nodes = topology.DescendantIdsOfNode( root_node_id );
      all_nodes.push_back( root_node_id );
      std::for_each( all_nodes.begin(), all_nodes.end(), [&topology, material]( std::uint64_t const id ) { topology.AddMaterialToNode( id, material ); } );
      topology.UpdateTopology();
   }
}

SCENARIO( "Topology offsets work properly on differently configured topologies", "[1rank]" ) {

   GIVEN( "A(ny) topology manager" ) {
      TopologyManager const topology = SingleRootNodeTopoloyWithMaximumLevel();
      WHEN( "The topology holds just one node" ) {
         REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );
         THEN( "The offset of rank zero is zero" ) {
            REQUIRE( topology.LeafOffsetOfRank( 0  ) == 0 );
         }
         THEN( "The offset of all other ranks is one" ) {
            REQUIRE( topology.LeafOffsetOfRank( 1  ) == 1 );
            REQUIRE( topology.LeafOffsetOfRank( 42 ) == 1 );
         }
      }
   }

   GIVEN( "A topology holding 1 parent and 8 children" ) {
      TopologyManager topology = SingleRootNodeTopoloyWithMaximumLevel( 1 );
      RefineZerothRootNode( topology );
      REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 9,8 ) );

      WHEN( "All nodes are on rank 0" ) {
         // That's the default.
         THEN( "The offset of rank zero is zero" ) {
            REQUIRE( topology.LeafOffsetOfRank( 0  ) == 0 );
         }
         THEN( "The offset of all other ranks is 8" ) {
            REQUIRE( topology.LeafOffsetOfRank( 1  ) == 8 );
            REQUIRE( topology.LeafOffsetOfRank( 42 ) == 8 );
         }
      }
      WHEN( "The nodes are distributed on 2 ranks" ) {
         // For load balancing to work, the nodes must have a weight = material inside them
         AddMaterialToAllNodes( topology, MaterialName::MaterialOne );
         topology.GetLoadBalancedTopology( 2 );

         THEN( "The offset of rank 0 is zero" ) {
            REQUIRE( topology.LeafOffsetOfRank( 0  ) == 0 );
         }
         THEN( "The offset of rank 1 is four" ) {
            REQUIRE( topology.LeafOffsetOfRank( 1  ) == 4 );
         }
         THEN( "The offset of the other ranks is eight" ) {
            REQUIRE( topology.LeafOffsetOfRank( 5  ) == 8 );
            REQUIRE( topology.LeafOffsetOfRank( 42 ) == 8 );
         }
      }
   }
}

SCENARIO( "The number of multi-phase nodes is correctly reported", "[1rank]" ) {
   GIVEN( "A topology with 8 Leaves on Lmax = 1 and one node on L0" ) {
      TopologyManager topology = SingleRootNodeTopoloyWithMaximumLevel( 2 );
      RefineZerothRootNode( topology );
      WHEN( "All nodes hold one material and the one leaf holds another one" ) {
         AddMaterialToAllNodes( topology, MaterialName::MaterialOne );
         topology.AddMaterialToNode( 0xA000000, MaterialName::MaterialTwo );
         topology.UpdateTopology();
         THEN( "The multi-node count is 1" ) {
            REQUIRE( topology.MultiPhaseNodeCount() == 1 );
         }
      }
      WHEN( "All nodes hold one material and 2 leaves and the root node hold another one" ) {
         AddMaterialToAllNodes( topology, MaterialName::MaterialOne );
         topology.AddMaterialToNode( root_node_id, MaterialName::MaterialTwo );
         topology.AddMaterialToNode( 0xA000001, MaterialName::MaterialTwo );
         topology.AddMaterialToNode( 0xA000002, MaterialName::MaterialTwo );
         topology.UpdateTopology();
         THEN( "The multi-node count is 3" ) {
            REQUIRE( topology.MultiPhaseNodeCount() == 3 );
         }
      }
   }
}