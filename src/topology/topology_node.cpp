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
#include "topology_node.h"

#include <algorithm>
#include "topology/id_information.h"
#include "user_specifications/compile_time_constants.h"
#include "utilities/container_operations.h"

namespace TNC = TopologyNodeConstants;

/**
 * @brief Constructs a topology node as leaf without children, without assigning a future rank and without materials.
 * @param rank Rank holding the node to be created.
 */
TopologyNode::TopologyNode( int const rank ) : current_rank_( rank ), target_rank_( TNC::unassigned_rank ), is_leaf_( true ), materials_() {}

/**
 * @brief Constructs a topology node as leaf with the given materials, but without assigning a future rank.
 * @param materials The materials to be present in the node.
 * @param rank Rank holding the node to be created.
 */
TopologyNode::TopologyNode( std::vector<MaterialName> const& materials, int const rank ) : current_rank_( rank ), target_rank_( TNC::unassigned_rank ), is_leaf_( true ), materials_( ContainerOperations::SortedCopy( materials ) ) {}

/**
 * @brief Adds the given material to the node.
 * @param material The material to be added to the node.
 * @note Does not perform checks if the material is already present. Inserting the same material twice will lead to undefined behavior.
 */
void TopologyNode::AddMaterial( MaterialName const material ) {
   materials_.insert( std::upper_bound( std::begin( materials_ ), std::end( materials_ ), material ), material );
}

/**
 * @brief Removes the given material from the node.
 * @param material The material to be removed from the node.
 */
void TopologyNode::RemoveMaterial( MaterialName const material ) {
   materials_.erase( std::remove( std::begin( materials_ ), std::end( materials_ ), material ), std::end( materials_ ) );
}

/**
 * @brief Gives the  materials present in the node.
 * @return The materials.
 */
std::vector<MaterialName> TopologyNode::Materials() const {
   return materials_;
}

/**
 * @brief Gives the material of a single phase node.
 * @return The material.
 * @note If called on a multi-phase node result is undefined.
 */
MaterialName TopologyNode::SingleMaterial() const {
   return materials_.front();
}

/**
 * @brief Gives the number of materials present in the node.
 * @return The number of materials.
 */
std::size_t TopologyNode::NumberOfMaterials() const {
   return materials_.size();
}

/**
 * @brief Converts a node into a parent.
 */
void TopologyNode::MakeParent() {
   is_leaf_ = false;
}

/**
 * @brief Converts a parent node into leaf.
 */
void TopologyNode::MakeLeaf() {
   is_leaf_ = true;
}

/**
 * @brief Whether the node is a leaf.
 * @return True if the node is a leaf, false otherwise.
 */
bool TopologyNode::IsLeaf() const {
   return is_leaf_;
}

/**
 * @brief Gives the rank the node resides on.
 * @return The rank id.
 */
int TopologyNode::Rank() const {
   return current_rank_;
}

/**
 * @brief Indicates if the node resides on the provided rank.
 * @return True if node is on the rank. False otherwise.
 */
bool TopologyNode::IsOnRank( int const rank ) const {
   return current_rank_ == rank;
}

/**
 * @brief Gives the rank on which the node should resides on (in the future).
 * @return The rank id.
 */
int TopologyNode::TargetRank() const {
   return target_rank_;
}

/**
 * @brief Sets the desired rank of the node to the given rank.
 * @param rank The rank the node should reside on.
 */
void TopologyNode::AssignTargetRank( int const rank ) {
   target_rank_ = rank;
}

/**
 * @brief Sets the rank of the node to the desired rank.
 */
void TopologyNode::SetCurrentRankAccordingToTargetRank() {
   current_rank_ = target_rank_;
}

/**
 * @brief Indicates whether the rank of the node is the desired rank.
 * @return True if the node's rank is the desired one. False otherwise.
 */
bool TopologyNode::IsBalanced() const {
   return current_rank_ == target_rank_;
}
