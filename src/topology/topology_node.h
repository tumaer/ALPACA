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
#ifndef TOPOLOGY_NODE_H
#define TOPOLOGY_NODE_H

#include <vector>
#include <tuple>
#include "topology/node_id_type.h"
#include "materials/material_definitions.h"

namespace TopologyNodeConstants {
   constexpr int unassigned_rank = -1;
}// namespace TopologyNodeConstants

/**
 * @brief The TopologyNode class organizes the light weight global ( over MPI ranks ) node information in a tree structure. Allowing the TopologyManager efficient searches.
 */
class TopologyNode {

   int current_rank_;
   int target_rank_;
   bool is_leaf_;
   std::vector<MaterialName> materials_;

public:
   explicit TopologyNode( int const rank = TopologyNodeConstants::unassigned_rank );
   explicit TopologyNode( std::vector<MaterialName> const& material, int const rank = TopologyNodeConstants::unassigned_rank );
   ~TopologyNode()                     = default;
   TopologyNode( TopologyNode const& ) = delete;
   TopologyNode& operator=( TopologyNode const& ) = delete;
   TopologyNode( TopologyNode&& )                 = delete;
   TopologyNode& operator=( TopologyNode&& ) = delete;

   void AddMaterial( MaterialName const material );
   void RemoveMaterial( MaterialName const material );

   std::vector<MaterialName> Materials() const;
   MaterialName SingleMaterial() const;
   std::size_t NumberOfMaterials() const;

   void MakeParent();
   void MakeLeaf();

   bool IsLeaf() const;

   int Rank() const;
   bool IsOnRank( int const rank ) const;
   int TargetRank() const;
   void AssignTargetRank( int const rank );
   void SetCurrentRankAccordingToTargetRank();

   bool IsBalanced() const;
};

#endif// TOPOLOGY_NODE_H
