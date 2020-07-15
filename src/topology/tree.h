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
#ifndef TREE_H
#define TREE_H

#include <vector>
#include <memory>
#include <unordered_map>
#include "topology_manager.h"
#include "block_definitions/interface_block.h"
#include "node.h"

/**
 * @brief The Tree class holds the information about local material data. Data on the MR levels is stored in Node containers as Binary-, Quad-, Oct-tree
 *        in one, two and three dimensions, respectively. Tree does not manipulate the data itself, it is an accessor that provides data to a
 *        solver or similar. No real tree-searches are performed due to the unique node indexing. Tree is rather a guarantee to access the correct
 *        data and that changes are transferable e.g. to the TopologyManager in the proper way. The tree must not change global data at any time.
 */
class Tree {
   // topology for global information of the nodes on all ranks
   TopologyManager const& topology_;
   // size of a cubic block on level zero (non-dimensionalized)
   double const node_size_on_level_zero_;
   // all nodes contained in this tree ( current rank )
   std::vector<std::unordered_map<nid_t, Node>> nodes_;

   void InsertNode( nid_t const id, std::vector<MaterialName> const materials, std::int8_t const interface_tag );

public:
   Tree() = delete;
   explicit Tree( TopologyManager const& topology, unsigned int const maximum_level, double const node_size_on_level_zero );
   ~Tree()             = default;
   Tree( Tree const& ) = delete;
   Tree& operator=( Tree const& ) = delete;
   Tree( Tree&& )                 = delete;
   Tree& operator=( Tree&& ) = delete;

   /**
    * @brief Gives the non-dimensionalized node size on level zero
    * @return Number of nodes
    */
   inline double GetNodeSizeOnLevelZero() const {
      return node_size_on_level_zero_;
   }

   // Functions for manipulating the tree
   Node& CreateNode( nid_t const id, std::vector<MaterialName> const& materials, std::int8_t const ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()],
                     std::unique_ptr<InterfaceBlock> interface_block = nullptr );
   Node& CreateNode( nid_t const id, std::vector<MaterialName> const& materials );
   void RemoveNodeWithId( nid_t const id );
   std::vector<nid_t> RefineNode( nid_t const id );

   // Functions to return leaf nodes
   std::vector<std::reference_wrapper<Node>> Leaves();
   std::vector<std::reference_wrapper<Node const>> Leaves() const;
   std::vector<std::reference_wrapper<Node>> LeavesOnLevel( unsigned int const level );
   std::vector<std::reference_wrapper<Node const>> LeavesOnLevel( unsigned int const level ) const;
   std::vector<std::reference_wrapper<Node>> NonLevelsetLeaves( unsigned int const level );
   std::vector<std::reference_wrapper<Node const>> NonLevelsetLeaves( unsigned int const level ) const;
   std::vector<std::reference_wrapper<Node const>> InterfaceLeaves() const;

   // Functions to return all nodes ( not only leaf nodes )
   std::vector<std::reference_wrapper<Node>> NodesOnLevel( unsigned int const level );
   std::vector<std::reference_wrapper<Node const>> NodesOnLevel( unsigned int const level ) const;
   std::vector<std::reference_wrapper<Node>> NodesWithLevelset();
   std::vector<std::reference_wrapper<Node const>> NodesWithLevelset() const;
   std::vector<std::reference_wrapper<Node const>> AllNodes() const;

   // Functions to return specific ndes
   Node const& GetNodeWithId( nid_t const id ) const;
   Node& GetNodeWithId( nid_t const id );
   std::pair<nid_t const, Node>& NodeIdPair( nid_t const id );
   std::pair<nid_t const, Node> const& NodeIdPair( nid_t const id ) const;
   std::unordered_map<nid_t, Node>& GetLevelContent( unsigned int const level );
   std::unordered_map<nid_t, Node> const& GetLevelContent( unsigned int const level ) const;

   /**
    * @brief Gives a reference to the complete node list in this tree instance, i e. the complete tree on current MPI rank.
    * @return List of nodes. An array for each level holding arbitrary number of nodes on each level.
    */
   inline std::vector<std::unordered_map<nid_t, Node>>& FullNodeList() { return nodes_; }
   /**
    * @brief Const overload.
    */
   inline std::vector<std::unordered_map<nid_t, Node>> const& FullNodeList() const { return nodes_; }
};

#endif// TREE_H
