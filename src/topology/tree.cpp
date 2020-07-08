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
#include "tree.h"

#include <stdexcept>
#include "topology/id_information.h"

/**
 * @brief Default constructor. Creates all nodes on level zero which should reside on the mpi rank the tree is on itself.
 * @param topology Reference to obtain information about the global state of the simulation.
 * @param maximum_level The maximum level possibly present in the simulation.
 * @param node_size_on_level_zero The geometric size of blocks on level zero.
 */
Tree::Tree( TopologyManager const& topology, unsigned int const maximum_level, double const node_size_on_level_zero ) : topology_( topology ),
                                                                                                                        node_size_on_level_zero_( node_size_on_level_zero ),
                                                                                                                        nodes_( maximum_level + 1 )// Level 0 + #Levels
{
   /** Empty besides initializer list */
}

/**
 * @brief Inserts a Node with specified id and materials into the tree structure.
 * @param id The unique id of the node to be created.
 * @param materials Identifiers of the materials present in the node to be created.
 * @param interface_tag The interface tag present in all of the node.
 */
void Tree::InsertNode( std::uint64_t const id, std::vector<MaterialName> const materials, std::int8_t const interface_tag ) {
   nodes_[LevelOfNode( id )].emplace( std::piecewise_construct, std::forward_as_tuple( id ), std::forward_as_tuple( id, node_size_on_level_zero_, materials, interface_tag ) );
}

/**
 * @brief Allows creation of nodes in this instance from the outside. Should therefore be used with extreme caution.
 * Needed e.g. at initalization.
 * @param id Unique identifier of the node to be created.
 * @param materials The materials to be contained in the new node.
 * @param interface_tags The interface tags to be contained in the new node.
 * @param interface_block The InterfaceBlock to be contained in the new node.
 */
Node& Tree::CreateNode( std::uint64_t const id, std::vector<MaterialName> const& materials,
                        std::int8_t const ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()], std::unique_ptr<InterfaceBlock> interface_block ) {

   unsigned int const level = LevelOfNode( id );
   auto entry_and_decision  = nodes_[level].emplace( std::piecewise_construct, std::forward_as_tuple( id ),
                                                    std::forward_as_tuple( id, node_size_on_level_zero_, materials, interface_tags, std::move( interface_block ) ) );

#ifndef PERFORMANCE
   if( !std::get<1>( entry_and_decision ) ) {
      throw std::logic_error( "Could not insert node into tree. Id already existed" );
   }
#endif
   return std::get<1>( *std::get<0>( entry_and_decision ) );
}

/**
 * @brief Overload, see also other implementation. Creates a node with default arguments for node constructor.
 * @param id .
 * @param materials .
 */
Node& Tree::CreateNode( std::uint64_t const id, std::vector<MaterialName> const& materials ) {

   unsigned int const level = LevelOfNode( id );
   auto entry_and_decision  = nodes_[level].emplace( std::piecewise_construct, std::forward_as_tuple( id ),
                                                    std::forward_as_tuple( id, node_size_on_level_zero_, materials ) );

#ifndef PERFORMANCE
   if( !std::get<1>( entry_and_decision ) ) {
      throw std::logic_error( "Could not insert node into tree. Id already existed" );
   }
#endif
   return std::get<1>( *std::get<0>( entry_and_decision ) );
}

/**
 * @brief Refines given Node in a three dimensional simulation, inserts its eight children into this tree instance and returns
 *        the ids of the created children. $NOT SAFE: Wrong input gives corrupted data/tree structure. Must only be called on leaves$.
 * @param id Id of the leaf which is to be refined, i.e. becomes a parent node.
 * @return List of childrens' ids.
 */
std::vector<std::uint64_t> Tree::RefineNode( std::uint64_t const id ) {

   unsigned int const level = LevelOfNode( id );
#ifndef PERFORMANCE
   if( level >= CC::AMNL() ) {
      throw std::invalid_argument( "Nodes on this level cannot be refined further" );
   }
#endif

   std::vector<std::uint64_t> const children_ids = IdsOfChildren( id );//IdsOfChildren adjusts to 1D/2D.
   Node const& node                              = GetNodeWithId( id );

   for( std::int64_t const child_id : children_ids ) {
      // only single material nodes are supposed to be refined, hence the use of uniform interface tag is valid
      InsertNode( child_id, topology_.GetMaterialsOfNode( id ), node.GetUniformInterfaceTag() );
   }

   return children_ids;
}

/**
 * @brief Removes the node with the given id.
 * @param id The identifier of the node to be removed.
 * @note Does not check for correctness. Node must be present otherwise undefiend behavior or exceptions will hunt you.
 */
void Tree::RemoveNodeWithId( std::uint64_t const id ) {
   unsigned int level = LevelOfNode( id );
   nodes_[level].erase( id );
}

/**
 * @brief Returns a list of all leaf nodes on this rank. $List is in arbitrary order$.
 * @return List of pointers to the leaves in this tree instance.
 */
std::vector<std::reference_wrapper<Node>> Tree::Leaves() {

   std::vector<std::uint64_t> leaf_ids = topology_.LocalLeafIds();
   std::vector<std::reference_wrapper<Node>> leaves;
   leaves.reserve( leaf_ids.size() );

   for( auto const& id : leaf_ids ) {
      leaves.emplace_back( GetNodeWithId( id ) );//We add this leaf
   }
   return leaves;
}

/**
 * @brief Const overload.
 */
std::vector<std::reference_wrapper<Node const>> Tree::Leaves() const {

   std::vector<std::uint64_t> leaf_ids = topology_.LocalLeafIds();
   std::vector<std::reference_wrapper<Node const>> leaves;
   leaves.reserve( leaf_ids.size() );

   for( auto const& id : leaf_ids ) {
      leaves.emplace_back( GetNodeWithId( id ) );//We add this leaf
   }
   return leaves;
}

/**
 * @brief Gives a list of all leaves on the specified level. $List is in arbitrary order$.
 * @param level The level of interest.
 * @return List of leaves.
 */
std::vector<std::reference_wrapper<Node>> Tree::LeavesOnLevel( unsigned int const level ) {

   std::vector<std::uint64_t> leaf_ids_on_level = topology_.LocalLeafIdsOnLevel( level );
   std::vector<std::reference_wrapper<Node>> leaves;
   leaves.reserve( leaf_ids_on_level.size() );

   for( auto& id : leaf_ids_on_level ) {
      leaves.emplace_back( GetNodeWithId( id ) );//We add this leaf
   }
   return leaves;
}

/**
* @brief const overload.
*/
std::vector<std::reference_wrapper<Node const>> Tree::LeavesOnLevel( unsigned int const level ) const {

   std::vector<std::uint64_t> leaf_ids_on_level = topology_.LocalLeafIdsOnLevel( level );
   std::vector<std::reference_wrapper<Node const>> leaves;
   leaves.reserve( leaf_ids_on_level.size() );

   for( auto const& id : leaf_ids_on_level ) {
      leaves.emplace_back( GetNodeWithId( id ) );//We add this leaf
   }
   return leaves;
}

/**
 * @brief Gives a list of leaves on the respective level which do not hold a levelset.
 * @param level The level of interest.
 * @return The requested node-list without levelset nodes
 */
std::vector<std::reference_wrapper<Node>> Tree::NonLevelsetLeaves( unsigned int const level ) {
   std::vector<std::uint64_t> const leaf_ids_on_level = topology_.LocalLeafIdsOnLevel( level );
   std::vector<std::reference_wrapper<Node>> non_levelset_leaves;
   non_levelset_leaves.reserve( leaf_ids_on_level.size() );

   for( auto const& id : leaf_ids_on_level ) {
      Node& node = GetNodeWithId( id );
      if( !node.HasLevelset() ) {
         non_levelset_leaves.emplace_back( node );
      }
   }

   return non_levelset_leaves;
}

/**
 * @brief Const overload.
 */
std::vector<std::reference_wrapper<Node const>> Tree::NonLevelsetLeaves( unsigned int const level ) const {
   std::vector<std::uint64_t> const leaf_ids_on_level = topology_.LocalLeafIdsOnLevel( level );
   std::vector<std::reference_wrapper<Node const>> non_levelset_leaves;
   non_levelset_leaves.reserve( leaf_ids_on_level.size() );

   for( auto const& id : leaf_ids_on_level ) {
      Node const& node = GetNodeWithId( id );
      if( !node.HasLevelset() ) {
         non_levelset_leaves.emplace_back( node );
      }
   }

   return non_levelset_leaves;
}

/**
 * @brief Gives a list of nodes that contain a levelset
 * @param level .
 * @return The requested node-list.
 */
std::vector<std::reference_wrapper<Node const>> Tree::InterfaceLeaves() const {

   std::vector<std::uint64_t> interface_leaf_ids = topology_.LocalInterfaceLeafIds();
   std::vector<std::reference_wrapper<Node const>> interface_leaves;
   interface_leaves.reserve( interface_leaf_ids.size() );

   for( auto const& id : interface_leaf_ids ) {
      interface_leaves.emplace_back( GetNodeWithId( id ) );//We add this leaf
   }
   return interface_leaves;
}

/**
 * @brief Gives a list of all nodes on the specified level. $List is in arbitrary order$.
 * @param level The level of interest.
 * @return List of nodes.
 */
std::vector<std::reference_wrapper<Node>> Tree::NodesOnLevel( unsigned int const level ) {

#ifndef PERFORMANCE
   if( level > nodes_.size() ) {
      throw std::invalid_argument( "Requested Level does not exist" );
   }
#endif

   std::vector<std::reference_wrapper<Node>> nodes;

   for( auto& node : nodes_[level] ) {
      nodes.emplace_back( node.second );
   }

   return nodes;
}

/**
 * @brief const overload.
 */
std::vector<std::reference_wrapper<Node const>> Tree::NodesOnLevel( unsigned int const level ) const {

#ifndef PERFORMANCE
   if( level > nodes_.size() ) {
      throw std::invalid_argument( "Requested Level does not exist" );
   }
#endif

   std::vector<std::reference_wrapper<Node const>> nodes;

   for( auto const& node : nodes_[level] ) {
      nodes.emplace_back( node.second );
   }

   return nodes;
}

/**
 * @brief gives all nodes with levelset in the tree. $List is in arbitrary order$.
 * @return List of nodes with levelset. List may be empty.
 */
std::vector<std::reference_wrapper<Node>> Tree::NodesWithLevelset() {
   //Levelset only on maximum level
   std::unordered_map<std::uint64_t, Node>& maximum_level = nodes_.back();
   std::vector<std::reference_wrapper<Node>> nodes;
   nodes.reserve( maximum_level.size() );

   for( auto& id_node : maximum_level ) {
      if( id_node.second.HasLevelset() ) {
         nodes.emplace_back( id_node.second );
      }
   }
   return nodes;
}

/**
 * @brief const overload.
 */
std::vector<std::reference_wrapper<Node const>> Tree::NodesWithLevelset() const {
   //Levelset only on maximum level
   std::unordered_map<std::uint64_t, Node> const& maximum_level = nodes_.back();
   std::vector<std::reference_wrapper<Node const>> nodes;
   nodes.reserve( maximum_level.size() );

   for( auto const& node : maximum_level ) {
      if( node.second.HasLevelset() ) {
         nodes.emplace_back( node.second );
      }
   }
   return nodes;
}

/**
 * @brief Gives all nodes in the tree
 * @return List of nodes.
 */
std::vector<std::reference_wrapper<Node const>> Tree::AllNodes() const {

   std::vector<std::uint64_t> node_ids = topology_.LocalNodeIds();
   std::vector<std::reference_wrapper<Node const>> nodes;
   nodes.reserve( node_ids.size() );

   for( auto const& id : node_ids ) {
      nodes.emplace_back( GetNodeWithId( id ) );//We add this node
   }
   return nodes;
}

/**
 * @brief Returns pointer to the node having the requested id. Throws exception if node is not in this tree,
 *        i.e. must only be called on correct mpi rank.
 * @param id Id to uniquely identify the Node.
 * @return Pointer to requested Node if existent. Throws exception otherwise.
 * @note Hotpath function.
 */
Node& Tree::GetNodeWithId( std::uint64_t const id ) {
   unsigned int const level = LevelOfNode( id );
   return nodes_[level].at( id );
}

/**
 * @brief const overload.
 * @note Hotpath function.
 */
Node const& Tree::GetNodeWithId( std::uint64_t const id ) const {
   unsigned int const level = LevelOfNode( id );
   return nodes_[level].at( id );
}

/**
 * @brief Gives access to entry with the requested id.
 * @param id .
 * @return .
 * @note Leads to undefined behaviour if called wiht non-exsiting id.
 */

std::pair<std::uint64_t const, Node>& Tree::NodeIdPair( std::uint64_t const id ) {
   unsigned int const level = LevelOfNode( id );
   return *nodes_[level].find( id );
}

/**
 * @brief const overload.
 */
std::pair<std::uint64_t const, Node> const& Tree::NodeIdPair( std::uint64_t const id ) const {
   unsigned int const level = LevelOfNode( id );
   return *nodes_[level].find( id );
}

/**
 * @brief Gives a list of all nodes in this instance on the specified level. List is in arbitrary order.
 * @param level The level of interest.
 * @return List of nodes.
 */
std::unordered_map<std::uint64_t, Node>& Tree::GetLevelContent( unsigned int const level ) {
#ifndef PERFORMANCE
   if( level > nodes_.size() ) {
      throw std::invalid_argument( "Requested Level does not exist" );
   }
#endif
   return nodes_[level];
}

/**
 * @brief const overload.
 */
std::unordered_map<std::uint64_t, Node> const& Tree::GetLevelContent( unsigned int const level ) const {
#ifndef PERFORMANCE
   if( level > nodes_.size() ) {
      throw std::invalid_argument( "Requested Level does not exist" );
   }
#endif
   return nodes_[level];
}