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
#include <catch2/catch.hpp>
#include "materials/material_definitions.h"
#include "topology/topology_node.h"

SCENARIO( "Topology Nodes created with different parameters, materials are modified, made parents and ranks are updated ", "[1rank]" ) {
   constexpr int test_rank                       = 42;
   std::vector<MaterialName> const two_materials = { MaterialName::MaterialOne, MaterialName::MaterialTwo };
   GIVEN( "Three topology nodes are constructed. An empty one, a rank constructed node and one with two materials." ) {
      TopologyNode empty_constructed    = TopologyNode();
      TopologyNode rank_constructed     = TopologyNode( test_rank );
      TopologyNode material_constructed = TopologyNode( two_materials );
      THEN( "The nodes are leaves" ) {
         REQUIRE( empty_constructed.IsLeaf() );
         REQUIRE( rank_constructed.IsLeaf() );
         REQUIRE( material_constructed.IsLeaf() );
      }
      THEN( "The rank constructed node is already on the given rank, the others are unassigned" ) {
         REQUIRE( rank_constructed.Rank() == test_rank );
         REQUIRE( empty_constructed.Rank() == TopologyNodeConstants::unassigned_rank );
         REQUIRE( material_constructed.Rank() == TopologyNodeConstants::unassigned_rank );
      }
      THEN( "The material constructed node already holds two materials, the others hold none" ) {
         REQUIRE( material_constructed.NumberOfMaterials() == two_materials.size() );
         REQUIRE( rank_constructed.NumberOfMaterials() == 0 );
         REQUIRE( empty_constructed.NumberOfMaterials() == 0 );
      }
      THEN( "Only the rank constructed yields true when asked if node is on the test rank" ) {
         REQUIRE( rank_constructed.IsOnRank( test_rank ) );
         REQUIRE_FALSE( empty_constructed.IsOnRank( test_rank ) );
         REQUIRE_FALSE( material_constructed.IsOnRank( test_rank ) );
      }

      WHEN( "We make the nodes parents followed by making them leaves" ) {
         empty_constructed.MakeParent();
         rank_constructed.MakeParent();
         material_constructed.MakeParent();
         THEN( "The nodes become non-leaves" ) {
            REQUIRE_FALSE( empty_constructed.IsLeaf() );
            REQUIRE_FALSE( rank_constructed.IsLeaf() );
            REQUIRE_FALSE( material_constructed.IsLeaf() );
         }
         empty_constructed.MakeLeaf();
         rank_constructed.MakeLeaf();
         material_constructed.MakeLeaf();
         THEN( "The nodes become leaves again" ) {
            REQUIRE( empty_constructed.IsLeaf() );
            REQUIRE( rank_constructed.IsLeaf() );
            REQUIRE( material_constructed.IsLeaf() );
         }
      }

      WHEN( "We add the first of the materials to the empty constructed node and remove the second material form the material constructed node" ) {
         empty_constructed.AddMaterial( two_materials.front() );
         material_constructed.RemoveMaterial( two_materials.back() );
         THEN( "Both nodes hold the same materials" ) {
            REQUIRE( empty_constructed.Materials() == material_constructed.Materials() );
         }
         THEN( "The single material on both is the first of the two materials" ) {
            REQUIRE( empty_constructed.SingleMaterial() == two_materials.front() );
            REQUIRE( material_constructed.SingleMaterial() == two_materials.front() );
         }
      }

      WHEN( "We assign the test rank as target rank to all nodes" ) {
         empty_constructed.AssignTargetRank( test_rank );
         rank_constructed.AssignTargetRank( test_rank );
         material_constructed.AssignTargetRank( test_rank );
         THEN( "This test_rank is correctly reported" ) {
            REQUIRE( empty_constructed.TargetRank() == test_rank );
            REQUIRE( rank_constructed.TargetRank() == test_rank );
            REQUIRE( material_constructed.TargetRank() == test_rank );
         }
         THEN( "Only the rank constructed node is balanced" ) {
            REQUIRE( rank_constructed.IsBalanced() );
            REQUIRE_FALSE( empty_constructed.IsBalanced() );
            REQUIRE_FALSE( material_constructed.IsBalanced() );
         }
      }

      WHEN( "We assigning a different rank to all nodes and update the ranks" ) {
         constexpr int other_rank = 2;
         empty_constructed.AssignTargetRank( other_rank );
         rank_constructed.AssignTargetRank( other_rank );
         material_constructed.AssignTargetRank( other_rank );
         empty_constructed.SetCurrentRankAccordingToTargetRank();
         rank_constructed.SetCurrentRankAccordingToTargetRank();
         material_constructed.SetCurrentRankAccordingToTargetRank();
         THEN( "The rank is reported as other rank" ) {
            REQUIRE( empty_constructed.Rank() == other_rank );
            REQUIRE( rank_constructed.Rank() == other_rank );
            REQUIRE( material_constructed.Rank() == other_rank );
         }
      }
   }
}
