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

#include "topology/tree.h"
#include "topology/topology_manager.h"

SCENARIO( "Retriving nodes from a tree with a different number of nodes on two levels", "[1rank]" ) {
   GIVEN( "A topology, a geometric size and a maximum level" ) {
      TopologyManager const topology       = TopologyManager();
      constexpr double geometric_size      = 1.0;
      constexpr unsigned int maximum_level = 1;
      WHEN( "An empty tree is created with these inputs" ) {
         Tree const tree = Tree( topology, geometric_size, maximum_level );
         THEN( "Retriving a node throws an exception" ) {
            REQUIRE_THROWS( tree.GetNodeWithId( 0 ) );
            REQUIRE_THROWS( tree.GetNodeWithId( topology.LeafIds().front() ) );
         }
      }
      WHEN( "When the tree holds one parent node" ) {
         Tree tree = Tree( topology, geometric_size, maximum_level );
         tree.CreateNode( 0x1400000, { MaterialName::MaterialOne } );
         THEN( "The number of nodes on level zero should contain one element and on level 1 zero elements" ) {
            REQUIRE( tree.NodesOnLevel( 0 ).size() == 1 );
            REQUIRE( tree.NodesOnLevel( 1 ).size() == 0 );
         }
      }
      WHEN( "When the tree holds one parent node and two child nodes" ) {
         Tree tree = Tree( topology, geometric_size, maximum_level );
         tree.CreateNode( 0x1400000, { MaterialName::MaterialOne } );
         tree.CreateNode( 0xA000000, { MaterialName::MaterialOne } );
         tree.CreateNode( 0xA000001, { MaterialName::MaterialOne } );
         THEN( "The number of nodes on level zero should contain one element and on level 1 two elements" ) {
            REQUIRE( tree.NodesOnLevel( 0 ).size() == 1 );
            REQUIRE( tree.NodesOnLevel( 1 ).size() == 2 );
         }
      }
   }
}