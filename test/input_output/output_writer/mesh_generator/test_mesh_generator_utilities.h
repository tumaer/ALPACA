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
#include "communication/mpi_utilities.h"

#include "topology/topology_manager.h"
#include "topology/tree.h"

namespace TestUtilities {
   /**
    * @brief Refines the first node in the topology
    * @param topology Topology that should be updated (indirect return)
    * @param do_load_balance Flag whether load balancign should be carried out
    */
   inline void RefineFirstNodeInTopology( TopologyManager& topology, bool const do_load_balance = false ) {
      REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );
      topology.RefineNodeWithId( 0x1400000 );
      topology.UpdateTopology();
      // For load balancing to work, the nodes must have a weight = materials inside them
      topology.AddMaterialToNode( 0x1400000, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000000, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000001, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000002, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000003, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000004, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000005, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000006, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000007, MaterialName::MaterialOne );
      topology.UpdateTopology();
      if( do_load_balance ) {
         REQUIRE( topology.GetMaximumLevel() >= 1 );
         topology.PrepareLoadBalancedTopology( MpiUtilities::NumberOfRanks() );
      }
   }

   /**
    * @brief Refines the two nodes in the topology
    * @param topology Topology that should be updated (indirect return)
    * @param do_load_balance Flag whether load balancign should be carried out
    */
   inline void RefineTwoNodesInTopology( TopologyManager& topology, bool const do_load_balance = false ) {
      REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );
      topology.RefineNodeWithId( 0x1400000 );
      topology.UpdateTopology();
      // For load balancing to work, the nodes must have a weight = materials inside them
      topology.AddMaterialToNode( 0x1400000, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000000, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000001, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000002, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000003, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000004, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000005, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000006, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0xA000007, MaterialName::MaterialOne );
      topology.UpdateTopology();
      // Refine one of the second nodes
      topology.RefineNodeWithId( 0xA000000 );
      topology.UpdateTopology();
      // For load balancing to work, the nodes must have a weight = materials inside them
      topology.AddMaterialToNode( 0x50000000, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0x50000001, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0x50000002, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0x50000003, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0x50000004, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0x50000005, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0x50000006, MaterialName::MaterialOne );
      topology.AddMaterialToNode( 0x50000007, MaterialName::MaterialOne );
      topology.UpdateTopology();

      if( do_load_balance ) {
         REQUIRE( topology.GetMaximumLevel() >= 2 );
         topology.PrepareLoadBalancedTopology( MpiUtilities::NumberOfRanks() );
      }
   }

   inline void GetFirstNodeInTopologyWithoutMultiLeaves( TopologyManager& topology, Tree& tree ) {
      REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );
      topology.UpdateTopology();
      topology.AddMaterialToNode( 0x1400000, MaterialName::MaterialOne );
      tree.CreateNode( 0x1400000, { MaterialName::MaterialOne } );
      topology.UpdateTopology();
   }

   inline void RefineFirstNodeInTopologyAndAddTwoMultiLeaves( TopologyManager& topology, Tree& tree ) {
      REQUIRE( topology.NodeAndLeafCount() == std::pair<unsigned int, unsigned int>( 1, 1 ) );
      topology.RefineNodeWithId( 0x1400000 );
      topology.UpdateTopology();

      unsigned int counter = 0;
      for( auto const id : topology.LocalLeafIds() ) {
         if( counter >= 2 ) {
            topology.AddMaterialToNode( id, MaterialName::MaterialOne );
            tree.CreateNode( id, { MaterialName::MaterialOne } );
         } else {
            topology.AddMaterialToNode( id, MaterialName::MaterialOne );
            topology.AddMaterialToNode( id, MaterialName::MaterialTwo );
            tree.CreateNode( id, { MaterialName::MaterialOne, MaterialName::MaterialTwo } );
         }
         counter++;
      }
      topology.UpdateTopology();
   }
}// namespace TestUtilities