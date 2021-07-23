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
#include "topology/space_filling_curve_order.h"
#include <vector>
#include "topology/id_information.h"
#include "topology/node_id_type.h"
#include "utilities/container_operations.h"
#include "topology/space_filling_curve_index.h"

namespace {
   /**
    * @brief Gives the ids of some arbitrary nodes which reside on the same level.
    */
   std::vector<nid_t> BunchOfNodesOnLevelThree() {
      auto const level_three_start = IdsOfChildren( IdsOfChildren( IdsOfChildren( IdSeed() ).front() ).front() ).front();
      return { level_three_start,
               EastNeighborOfNodeWithId( level_three_start ),
               NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( level_three_start ) ),
               TopNeighborOfNodeWithId( NorthNeighborOfNodeWithId( EastNeighborOfNodeWithId( level_three_start ) ) ),
               TopNeighborOfNodeWithId( NorthNeighborOfNodeWithId( TopNeighborOfNodeWithId( level_three_start ) ) ),
               NorthNeighborOfNodeWithId( TopNeighborOfNodeWithId( level_three_start ) ),
               TopNeighborOfNodeWithId( level_three_start ) };
   }
}// namespace

SCENARIO( "Space-filling curves give the correct ordering of ids", "[1rank]" ) {
   GIVEN( "A bunch of ids on the same level" ) {
      auto nodes = BunchOfNodesOnLevelThree();
      WHEN( "We ask for the Hilbert Curve order of these ids and on a rotated input" ) {
         auto rotated_nodes = ContainerOperations::RotatedLeftCopy( nodes, 4 );
         OrderNodeIdsBySpaceFillingCurve( nodes, HilbertIndex );
         OrderNodeIdsBySpaceFillingCurve( rotated_nodes, HilbertIndex );
         THEN( "The mappings are identical" ) {
            REQUIRE( nodes == rotated_nodes );
         }
         THEN( "The mappings are unique" ) {
            REQUIRE( ContainerOperations::HoldsUniqueElements( nodes ) );
         }
      }
      WHEN( "We ask for the Z-Curve order these ids and on a rotated input" ) {
         auto rotated_nodes = ContainerOperations::RotatedLeftCopy( nodes, 4 );
         OrderNodeIdsBySpaceFillingCurve( nodes, LebesgueIndex );
         OrderNodeIdsBySpaceFillingCurve( rotated_nodes, LebesgueIndex );
         THEN( "The mappings are identical" ) {
            REQUIRE( nodes == rotated_nodes );
         }
         THEN( "The mapping are unique" ) {
            REQUIRE( ContainerOperations::HoldsUniqueElements( rotated_nodes ) );
         }
      }
   }
}
