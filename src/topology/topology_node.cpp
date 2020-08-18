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
#include "topology_node.h"

#include <algorithm>
#include "topology/id_information.h"
#include "user_specifications/compile_time_constants.h"
#include "utilities/container_operations.h"

/**
 * @brief Standard constructor for a TopologyNode. Initializes the node as a leaf without children and without assigning a future rank.
 * @param id The id of the node to be created.
 * @param rank Rank holding the node to be created.
 */
TopologyNode::TopologyNode( nid_t const id, int const rank ) : unique_id_( id ),
                                                               current_rank_( rank ),
                                                               future_rank_( -1 ),
                                                               children_{},
                                                               is_leaf_( true ),
                                                               materials_{} {
   // Empty Constructor, besides initializer list.
}

/**
 * @brief Traces the child with the given id and triggers refinement.
 * @param id The id of the node that is to be refined.
 */
void TopologyNode::Refine( nid_t const id ) {
   if( unique_id_ == id ) {
      if( is_leaf_ ) {
         Refine();
      }
   } else {
      GetChildWithId( id ).Refine( id );
   }
}

/**
 * @brief Refines the current leaf, i.e. inserts children as leaves and changes status away from leaf.
 */
void TopologyNode::Refine() {
   is_leaf_         = false;
   nid_t working_id = unique_id_ << 3;
   for( unsigned int i = 0; i < CC::NOC(); i++ ) {
      children_.emplace_back( TopologyNode( working_id + i, current_rank_ ) );
   }
   // DO NOT add materials to child nodes as this will be done later and not all materials might be present in the child
}

/**
 * @brief Coarsens the node with the given id. I.e. removes all its descendants and sets its status to leaf.
 * @param id The id of the node to be made a leaf.
 */
void TopologyNode::Coarse( nid_t const id ) {
   if( unique_id_ == id ) {
      for( auto& child : children_ ) {
         child.Coarse( child.Id() );
      }
      is_leaf_ = true;
      children_.clear();
   } else {
      if( !is_leaf_ ) {//coarsening might happen on more than one level, do nothing if already coarsened
         GetChildWithId( id ).Coarse( id );
      }
   }
}

/**
 * @brief Gives the node and leaf count, i.e. how many nodes and leaves are descending from this one.
 * @param nodes_and_leaves Indirect return parameter of the final result.
 */
void TopologyNode::NodeLeafCount( std::pair<unsigned int, unsigned int>& nodes_and_leaves ) const {
   nodes_and_leaves.first += 1;
   if( is_leaf_ )
      nodes_and_leaves.second += 1;
   else {
      for( TopologyNode const& child : children_ ) {
         child.NodeLeafCount( nodes_and_leaves );
      }
   }
}

/**
 * @brief Gives the amount of interface leaves within the decendants (itself included) of this node.
 * @brief retrun The interface node count.
 */
unsigned int TopologyNode::InterfaceLeafCount() const {
   unsigned int count = 0;
   if( is_leaf_ ) {
      if( materials_.size() != 1 ) {
         count++;
      }
   } else {
      for( TopologyNode const& child : children_ ) {
         count += child.InterfaceLeafCount();
      }
   }
   return count;
}

/**
 * @brief Gives the count of nodes and leaves in a list where each element corresponds to one MPI rank
 * @param nodes_leaves_per_rank Indirect return parameter.
 */
void TopologyNode::RankWiseNodeLeafCount( std::vector<std::pair<unsigned int, unsigned int>>& nodes_leaves_per_rank ) const {

   if( static_cast<int>( nodes_leaves_per_rank.size() ) < current_rank_ + 1 ) {// length starts with 0, int due to MPI
      nodes_leaves_per_rank.resize( current_rank_ + 1, std::make_pair( 0, 0 ) );
   }

   nodes_leaves_per_rank[current_rank_].first += 1;
   if( is_leaf_ ) {
      nodes_leaves_per_rank[current_rank_].second += 1;
   } else {
      for( TopologyNode const& child : children_ ) {
         child.RankWiseNodeLeafCount( nodes_leaves_per_rank );
      }
   }
}

/**
 * @brief Gives the count of interface leaves in a list where each element corresponds to one MPI rank.
 * @param nodes_leaves_per_rank Indirect return parameter.
 */
void TopologyNode::RankWiseInterfaceLeafCount( std::vector<unsigned int>& nodes_leaves_per_rank ) const {

   if( static_cast<int>( nodes_leaves_per_rank.size() ) < current_rank_ + 1 ) {// length starts with 0, int due to MPI
      nodes_leaves_per_rank.resize( current_rank_ + 1, 0u );
   }

   if( is_leaf_ ) {
      if( materials_.size() != 1 ) {
         nodes_leaves_per_rank[current_rank_] += 1;
      }
   } else {
      for( TopologyNode const& child : children_ ) {
         child.RankWiseInterfaceLeafCount( nodes_leaves_per_rank );
      }
   }
}

/**
 * @brief Gives the node and block count, i.e. how many nodes and blocks are descending from this one.
 * @param nodes_and_blocks Indirect return parameter of the final result.
 */
void TopologyNode::NodeBlockCount( std::pair<unsigned int, unsigned int>& nodes_and_blocks ) const {
   nodes_and_blocks.first += 1;
   nodes_and_blocks.second += materials_.size();
   if( !is_leaf_ ) {
      for( TopologyNode const& child : children_ ) {
         child.NodeBlockCount( nodes_and_blocks );
      }
   }
}
/**
 * @brief Gives the number of multi-phase nodes among the descendants of this node.
 * @return Number of multiphase nodes among children
 */
unsigned int TopologyNode::MultiPhaseNodeCount() const {
   unsigned int count = 0;
   for( TopologyNode const& child : children_ ) {
      count += child.MultiPhaseNodeCount();
   }
   auto materials_in_node_count = materials_.size();
   if( materials_in_node_count > 1 ) {
      return count + 1;
   } else {
      return count;
   }
}

/**
 * @brief Gives the count of nodes and blocks in a list where each element corresponds to one MPI rank.
 * @param nodes_blocks_per_rank Indirect return parameter.
 */
void TopologyNode::RankWiseNodeBlockCount( std::vector<std::pair<unsigned int, unsigned int>>& nodes_blocks_per_rank ) const {

   if( static_cast<int>( nodes_blocks_per_rank.size() ) < current_rank_ + 1 ) {// length starts with 0, int due to MPI
      nodes_blocks_per_rank.resize( current_rank_ + 1, std::make_pair( 0, 0 ) );
   }

   nodes_blocks_per_rank[current_rank_].first += 1;
   nodes_blocks_per_rank[current_rank_].second += materials_.size();
   if( !is_leaf_ ) {
      for( TopologyNode const& child : children_ ) {
         child.RankWiseNodeBlockCount( nodes_blocks_per_rank );
      }
   }
}

/**
 * @brief Gives the depth of the children, i.e. the number of generations in the tree aka. the number of refinement levels.
 * @return Depth of node
 * @note The tree depth starts with 1 by definition.
 */
unsigned int TopologyNode::GetDepth() const {
   if( is_leaf_ ) {
      return 1;
   } else {
      std::vector<unsigned int> child_depth;
      child_depth.reserve( CC::NOC() );
      for( TopologyNode const& child : children_ ) {
         child_depth.emplace_back( child.GetDepth() );
      }
      return 1 + *std::max_element( child_depth.begin(), child_depth.end() );
   }
}

/**
 * @brief Collects the ids of all nodes on the given level.
 * @param level The level of interest
 * @param ids Indirect return parameter
 */
void TopologyNode::IdsOnLevel( unsigned int const level, std::vector<nid_t>& ids ) const {
   if( LevelOfNode( unique_id_ ) == level ) {
      ids.emplace_back( unique_id_ );
   } else {
      for( TopologyNode const& child : children_ ) {
         child.IdsOnLevel( level, ids );
      }
   }
}

/**
 * @brief Gives the ids of all nodes on the given level that are handled by the current rank.
 * @param level Level of interest.
 * @param ids Indirect return parameter.
 */
void TopologyNode::LocalIdsOnLevel( unsigned int const level, std::vector<nid_t>& ids, int const local_rank ) const {
   if( LevelOfNode( unique_id_ ) == level && current_rank_ == local_rank ) {
      ids.emplace_back( unique_id_ );
   } else {
      for( TopologyNode const& child : children_ ) {
         child.LocalIdsOnLevel( level, ids, local_rank );
      }
   }
}

/**
 * @brief Gives the rank of the node with the given id.
 * @param id The id of the node of interest.
 * @return Current rank of the node.
 */
int TopologyNode::GetRank( nid_t const id ) const {
   if( unique_id_ == id )
      return current_rank_;
   else
      return GetChildWithId( id ).GetRank( id );
}

/**
 * @brief Gives a list of all leaves in this Topology node on the rank of interest.
 * @param local_leaves Indirect return parameter: The found local leaves.
 * @param rank The id of the local rank.
 * @return The number of found local leaves.
 */
unsigned int TopologyNode::LocalLeaves( std::vector<nid_t>& local_leaves, int const rank ) const {
   if( is_leaf_ ) {
      if( current_rank_ == rank ) {
         local_leaves.push_back( unique_id_ );
         return 1;
      } else {
         return 0;
      }
   } else {
      int sum = 0;
      for( TopologyNode const& child : children_ ) {
         sum += child.LocalLeaves( local_leaves, rank );
      }
      return sum;
   }
}

/**
 * @brief Gives a list of all leaves containing a levelset on the rank of interest
 * @param local_leaves Indirect return parameter: The found local leaves.
 * @param rank The id of the local rank.
 * @return The number of found local leaves.
 */
unsigned int TopologyNode::LocalInterfaceLeaves( std::vector<nid_t>& local_leaves, int const rank ) const {
   if( is_leaf_ ) {
      if( current_rank_ == rank && materials_.size() != 1 ) {
         local_leaves.push_back( unique_id_ );
         return 1;
      } else {
         return 0;
      }
   } else {
      int sum = 0;
      for( TopologyNode const& child : children_ ) {
         sum += child.LocalInterfaceLeaves( local_leaves, rank );
      }
      return sum;
   }
}

/**
 * @brief Gives a list of all nodes in this Topology node on the rank of interest.
 * @param local_nodes Indirect return parameter: The found local nodes.
 * @param rank The id of the local rank.
 */
void TopologyNode::LocalNodes( std::vector<nid_t>& local_nodes, int const rank ) const {
   // Add always this node
   if( current_rank_ == rank ) {
      local_nodes.push_back( unique_id_ );
   }
   // Loop through all children and add them subsequently
   for( TopologyNode const& child : children_ ) {
      child.LocalNodes( local_nodes, rank );
   }
}

/**
 * @brief Gives a list of all leaves in this TopologyNode ( tree ).
 * @param leaves Indirect return parameter: The found leaves.
 * @return The number of leaves.
 */
unsigned int TopologyNode::GetLeafIds( std::vector<nid_t>& leaves ) const {
   if( is_leaf_ ) {
      leaves.push_back( unique_id_ );
      return 1;
   } else {
      int sum = 0;
      for( TopologyNode const& child : children_ ) {
         sum += child.GetLeafIds( leaves );
      }
      return sum;
   }
}

/**
 * @brief Gives a list of all leaves on the specified rank in this TopologyNode ( tree ).
 * @param local_leaves Indirect return parameter: The found local leaves.
 * @param rank The id of the local rank.
 * @param level The level of interest
 * @param currentLevel The recursion level.
 * @return The number of local leaves.
 */
unsigned int TopologyNode::LocalLeavesOnLevel( std::vector<nid_t>& local_leaves, int const rank, unsigned int const level, unsigned int const current_level ) const {
   if( current_level == level ) {
      if( is_leaf_ && current_rank_ == rank ) {
         local_leaves.push_back( unique_id_ );
         return 1;
      } else {
         return 0;
      }
   } else {
      int sum = 0;
      for( TopologyNode const& child : children_ ) {
         sum += child.LocalLeavesOnLevel( local_leaves, rank, level, current_level + 1 );
      }
      return sum;
   }
}

/**
 * @brief Gives a list of all leaves in this TopologyNode ( tree ).
 * @param leaves Indirect return parameter: The found leaves.
 * @param current_level The recursion level.
 * @return The number of leaves.
 */
unsigned int TopologyNode::GetLeafIdsOnLevel( std::vector<nid_t>& leaves, unsigned int const level, unsigned int const current_level ) const {
   if( current_level == level ) {
      if( is_leaf_ ) {
         leaves.push_back( unique_id_ );
         return 1;
      } else {
         return 0;
      }
   } else {
      int sum = 0;
      for( TopologyNode const& child : children_ ) {
         sum += child.GetLeafIdsOnLevel( leaves, level, current_level + 1 );
      }
      return sum;
   }
}

/**
 * @brief Checks whether or not the node with the given id exists in the tree.
 * @param id The id of interest.
 * @return True if the node exists, false otherwise ( If all children respond with false ).
 */
bool TopologyNode::NodeExists( nid_t const id ) const {
   if( unique_id_ == id )
      return true;
   else {
      if( is_leaf_ ) {
         return false;
      }
      return GetChildWithId( id ).NodeExists( id );
   }
}

/**
 * @brief Checks whether or not the node with the given id is a leaf in the tree.
 * @param id The id of interest.
 * @return True if the node is a leaf, false otherwise.
 */
bool TopologyNode::NodeIsLeaf( nid_t const id ) const {
   if( unique_id_ == id )
      return is_leaf_;
   else {
      return GetChildWithId( id ).NodeIsLeaf( id );
   }
}

/**
 * @brief Gives a reference to the node with the given id.
 * @param id The id of the desired node.
 * @return The node with given id.
 * @note Does not directly return the correct value, gives a value for recursive search. I.e. if a grand-child is searched this function returns the mother of the grandchild.
 */
TopologyNode& TopologyNode::GetChildWithId( nid_t const id ) {
   // We have to search through the ancestor of the searched node until we find the child of the current one (recursive calls, remember...)
   nid_t id_on_child_level = id;
   while( LevelOfNode( id_on_child_level ) > ( LevelOfNode( unique_id_ ) + 1 ) ) {
      id_on_child_level = ParentIdOfNode( id_on_child_level );
   }
   return children_[PositionOfNodeAmongSiblings( id_on_child_level )];
}

/**
 * @brief Const variant of GetChildWithId(nid_t const id), see there for details.
 */
TopologyNode const& TopologyNode::GetChildWithId( nid_t const id ) const {
   // We have to search through the ancestor of the searched node until we find the child of the current one (recursive calls, remember...)
   nid_t id_on_child_level = id;
   while( LevelOfNode( id_on_child_level ) > ( LevelOfNode( unique_id_ ) + 1 ) ) {
      id_on_child_level = ParentIdOfNode( id_on_child_level );
   }
   return children_[PositionOfNodeAmongSiblings( id_on_child_level )];
}

/**
 * @brief Ensures that the load ( sum of weights ) is evenly distributed between the ranks
 *        Assigns the balanced topology to the future_rank_ member, which is then to be used in the actual redistribution in the LoadBalancing routine.
 * @return The rank to which this node should be moved.
 * @note Returns the modified member as of usage in recursive call.
 */
int TopologyNode::BalanceTargetRanks() {

   if( is_leaf_ ) {
#ifndef PERFORMANCE
      if( future_rank_ >= 0 ) {
         return future_rank_;
      } else {
         throw std::logic_error( "TopologyNode::balance: Rank for leaves must be set!" );
      }
#else
      return future_rank_;
#endif

   } else {
      std::vector<int> child_ranks;
      child_ranks.reserve( CC::NOC() );
      for( TopologyNode& child : children_ ) {
         child_ranks.push_back( child.BalanceTargetRanks() );
      }
      future_rank_ = ContainerOperations::MostFrequentElement( child_ranks );
   }
   return future_rank_;
}

/*
 * @brief Assigns the target rank to the leaf with the given id.
 * @param id The lead id to be assigned a new rank.
 * @param rank The rank to hold the node in the future.
 * @note Does not change the current rank of the leaf.
 */
void TopologyNode::AssignTargetRankToLeaf( nid_t const id, int const rank ) {
   if( id == unique_id_ ) {
#ifndef PERFORMANCE
      if( !is_leaf_ ) {
         throw std::logic_error( "Avoid calling assigning ranks to non-leaves!" );
      }
#endif
      future_rank_ = rank;
   } else {
      GetChildWithId( id ).AssignTargetRankToLeaf( id, rank );
   }
}

/**
 * @brief Gives a list of all nodes which are to be moved from one MPI rank to another in order to achieve a better load balance.
 * @param ids_current_future_rank_map Indirect return parameter, giving a list of all nodes that need to be shifted from one mpi rank to another.
 *        The entries are as the name suggest id - current rank - future rank.
 */
void TopologyNode::ListUnbalancedNodes( std::vector<std::tuple<nid_t const, int const, int const>>& ids_current_future_rank_map ) {
   if( future_rank_ != current_rank_ ) {
      ids_current_future_rank_map.push_back( std::make_tuple( unique_id_, current_rank_, future_rank_ ) );
      current_rank_ = future_rank_;
   }

   if( !is_leaf_ ) {
      for( TopologyNode& child : children_ ) {
         child.ListUnbalancedNodes( ids_current_future_rank_map );
      }
   }
}

/**
 * @brief Traces the child with the given id and adds the given material.
 * @param id The id of the node the material should be added to.
 * @param material The material to be added to the node.
 */
void TopologyNode::AddMaterial( nid_t const id, MaterialName const material ) {
   if( unique_id_ == id ) {
      AddMaterial( material );
   } else {
      GetChildWithId( id ).AddMaterial( id, material );
   }
}

/**
 * @brief Adds the given material to the node.
 * @param material The material to be added to the node.
 */
void TopologyNode::AddMaterial( MaterialName const material ) {
   materials_.emplace_back( material );
   std::sort( materials_.begin(), materials_.end(), []( MaterialName const a, MaterialName const b ) { return MTI( a ) < MTI( b ); } );
}

/**
 * @brief Traces the child with the given id and removes the given material.
 * @param id The id of the node the material should be removed from.
 * @param material The material to be removed from the node.
 */
void TopologyNode::RemoveMaterial( nid_t const id, MaterialName const material ) {
   if( unique_id_ == id ) {
      RemoveMaterial( material );
   } else {
      GetChildWithId( id ).RemoveMaterial( id, material );
   }
}

/**
 * @brief Removes the given material from the node.
 * @param material The material to be removed from the node.
 */
void TopologyNode::RemoveMaterial( MaterialName const material ) {
   materials_.erase( std::remove( materials_.begin(), materials_.end(), material ), materials_.end() );
   std::sort( materials_.begin(), materials_.end(), []( MaterialName const a, MaterialName const b ) { return MTI( a ) < MTI( b ); } );
}

/**
 * @brief Traces the child with the given id and returns its materials.
 * @param id The id of the node.
 * @return The materials of the given node.
 */
std::vector<MaterialName> TopologyNode::GetMaterials( nid_t const id ) const {
   return unique_id_ == id ? GetMaterials() : GetChildWithId( id ).GetMaterials( id );
}

/**
 * @brief Returns the materials of this node.
 * @return The materials of this node.
 */
std::vector<MaterialName> TopologyNode::GetMaterials() const {
   return materials_;
}

/**
 * @brief Returns the single material of a node based on its unique id.
 * @return The material name.
 */
MaterialName TopologyNode::GetSingleMaterial( nid_t const id ) const {
   return unique_id_ == id ? GetSingleMaterial() : GetChildWithId( id ).GetSingleMaterial( id );
}

/**
 * @brief Gives the material of a single phase node.
 * @return The material.
 */
MaterialName TopologyNode::GetSingleMaterial() const {
#ifndef PERFORMANCE
   if( materials_.size() != 1 ) {
      throw std::logic_error( "Cannot find single material in a multiphase node" );
   }
#endif
   return materials_.front();
}