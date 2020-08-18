
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
#include "topology_manager.h"

#include <algorithm>
#include <utility>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <numeric>
#include <functional>
#include <vector>

#include "communication/mpi_utilities.h"
#include "topology/id_information.h"
#include "topology/node_id_type.h"
#include "topology/space_filling_curve_order.h"
#include "utilities/string_operations.h"

namespace {
   /**
    * @brief Gives a count of elements that should be on each respective rank to obtained a balanced load.
    * @param number_of_elements The amount of elements to distribute accross the ranks.
    * @param number_of_ranks The amount of ranks to distribute the elements onto.
    * @return Vector of size number_of_ranks. Each entry gives the amount of elements the respective rank should hold in a well-balanced scenario.
    */
   std::vector<std::size_t> ElementsPerRank( std::size_t const number_of_elements, int const number_of_ranks ) {
      // Cast is save. Negative ranks indicate something is broken already anyways.
      std::size_t rank_count = static_cast<std::size_t>( number_of_ranks );
      std::size_t remainder  = number_of_elements % rank_count;
      std::vector<std::size_t> elements_per_rank( rank_count );
      for( std::size_t i = 0; i < rank_count; ++i ) {
         std::size_t fraction = 0;
         fraction += number_of_elements / rank_count;
         if( i < remainder ) {
            fraction++;
         }
         elements_per_rank[i] = fraction;
      }
      return elements_per_rank;
   }
}// namespace

/**
 * @brief Default constructor. Creates the local and global Id list on level zero. Default arguments allow testability.
 * @param maximum_level The possible maximum level in this simulation.
 * @param level_zero_nodes_x, level_zero_nodes_y, level_zero_nodes_z Number of blocks on level zero in the x/y/z-axis extension.
 * @param active_periodic_locations Side of the domain on which periodic boundaries are activated.
 */
TopologyManager::TopologyManager( std::array<unsigned int, 3> const level_zero_blocks, unsigned int const maximum_level, unsigned int const active_periodic_locations ) : maximum_level_( maximum_level ),
                                                                                                                                                                          active_periodic_locations_( active_periodic_locations ),
                                                                                                                                                                          number_of_nodes_on_level_zero_( level_zero_blocks ),
                                                                                                                                                                          forest_{},
                                                                                                                                                                          coarsenings_since_load_balance_{ 0 },
                                                                                                                                                                          refinements_since_load_balance_{ 0 } {
   nid_t id = IdSeed();

   std::vector<nid_t> initialization_list;

   /*
    * The Nodes are created in a spatial fashion traversing through X, Y and finally Z.
    * The ids in this traversal are not continuous due to the implicit shadow levels in the tree
    * Therefore some non-straightforward index magic needs to be applied
    * Initialization happens only on Level 0
    */
   for( unsigned int i = 0; i < number_of_nodes_on_level_zero_[2]; ++i ) {
      for( unsigned int j = 0; j < number_of_nodes_on_level_zero_[1]; ++j ) {
         for( unsigned int k = 0; k < number_of_nodes_on_level_zero_[0]; ++k ) {
            forest_.emplace_back( id, 0 );
            initialization_list.push_back( id );
            id = EastNeighborOfNodeWithId( initialization_list.back() );//find eastern neighbor of the just created Node
         }
         // Index magic to create the correct node once the ( outer ) loop counters are resetted
         id = NorthNeighborOfNodeWithId( initialization_list[number_of_nodes_on_level_zero_[0] * ( j + number_of_nodes_on_level_zero_[1] * i )] );
      }
      id = TopNeighborOfNodeWithId( initialization_list[( number_of_nodes_on_level_zero_[1] * number_of_nodes_on_level_zero_[0] * i )] );
   }

   //Topology Tree Node creation
   forest_.shrink_to_fit();

   // Assign correct ranks to nodes
   PrepareLoadBalancedTopology( MpiUtilities::NumberOfRanks() );
}

/**
 * @brief Updates the topology based on the recorded refinements, coarsening, weights and material changes.
 * @return True if Communication_managers cache needs to be invalidated.
 * @note The coarsening list is specially guarded, it needs to be flushed before it is considered here.
 */
bool TopologyManager::UpdateTopology() {

   // Flag that specifies whether the communication manager should update its neighbor relations
   bool invalidate_communication_manager_cache = false;
   int const number_of_ranks                   = MpiUtilities::NumberOfRanks();

   //Tree update
   //refine
   std::vector<nid_t> global_refine_list;
   MpiUtilities::LocalToGlobalData( local_refine_list_, MPI_LONG_LONG_INT, number_of_ranks, global_refine_list );

   for( auto const& refine_id : global_refine_list ) {
      forest_[PositionOfNodeInZeroTopology( refine_id )].Refine( refine_id );
   }
   // Invalididate cache if any node has been refined
   if( global_refine_list.size() > 0 ) {
      invalidate_communication_manager_cache = true;
   }

   local_refine_list_.clear();

   refinements_since_load_balance_ += global_refine_list.size();

   //UPDATE MATERIALS OF NODES
   std::tuple<std::vector<nid_t>, std::vector<MaterialName>> global_materials_list;
   MpiUtilities::LocalToGlobalData( std::get<0>( local_added_materials_list_ ), MPI_LONG_LONG_INT, number_of_ranks, std::get<0>( global_materials_list ) );
   MpiUtilities::LocalToGlobalData( std::get<1>( local_added_materials_list_ ), MPI_UNSIGNED_SHORT, number_of_ranks, std::get<1>( global_materials_list ) );

#ifndef PERFORMANCE
   if( std::get<0>( global_materials_list ).size() != std::get<1>( global_materials_list ).size() ) {
      throw std::logic_error( "Thou shall not create material-add-lists of unequal length" );
   }
#endif

   for( unsigned int i = 0; i < std::get<0>( global_materials_list ).size(); ++i ) {
      int const position = PositionOfNodeInZeroTopology( std::get<0>( global_materials_list )[i] );

      if( position >= 0 && position < static_cast<int>( forest_.size() ) ) {
         // Forest size cannot exceed 128x128x128 = 10^6, fits in int ( 10^9 positive values ), cast is safe.
         forest_[position].AddMaterial( std::get<0>( global_materials_list )[i], std::get<1>( global_materials_list )[i] );
      }
#ifndef PERFORMANCE
      else {
         throw std::logic_error( "TopologyManager::UpdateIdsList root tree does not exist!" );
      }
#endif
   }

   std::get<0>( local_added_materials_list_ ).clear();
   std::get<1>( local_added_materials_list_ ).clear();

   // We reuse the global list
   std::get<0>( global_materials_list ).clear();
   std::get<1>( global_materials_list ).clear();

   MpiUtilities::LocalToGlobalData( std::get<0>( local_removed_materials_list_ ), MPI_LONG_LONG_INT, number_of_ranks, std::get<0>( global_materials_list ) );
   MpiUtilities::LocalToGlobalData( std::get<1>( local_removed_materials_list_ ), MPI_UNSIGNED_SHORT, number_of_ranks, std::get<1>( global_materials_list ) );

#ifndef PERFORMANCE
   if( std::get<0>( global_materials_list ).size() != std::get<1>( global_materials_list ).size() ) {
      throw std::logic_error( "Created list of unequal length" );
   }
#endif

   for( unsigned int i = 0; i < std::get<0>( global_materials_list ).size(); ++i ) {
      int const position = PositionOfNodeInZeroTopology( std::get<0>( global_materials_list )[i] );
      // Forest size cannot exceed 128x128x128 = 10^6, fits in int ( 10^9 positive values ), cast is safe.
      if( position >= 0 && position < static_cast<int>( forest_.size() ) ) {
         forest_[position].RemoveMaterial( std::get<0>( global_materials_list )[i], std::get<1>( global_materials_list )[i] );
      }
#ifndef PERFORMANCE
      else {
         throw std::logic_error( "TopologyManager::UpdateIdsList root tree does not exist!" );
      }
#endif
   }

   std::get<0>( local_removed_materials_list_ ).clear();
   std::get<1>( local_removed_materials_list_ ).clear();

   return invalidate_communication_manager_cache;
}

/**
 * @brief Marks the node with the given id for refinement.
 * @param id The id of the leaf that is to be refined.
 * @note The actual refinement happens in a bundled fashion in another function. $No checks for leaves are performed caller is responsible$
 */
void TopologyManager::RefineNodeWithId( nid_t const id ) {
   local_refine_list_.push_back( id );
}

/**
 * @brief Adds the parent whose children may be coarsened to the coarsening list. The actual data deletion is than bundled using this list.
 * @param parent_id The id of the node that is to be made a leaf.
 */
void TopologyManager::CoarseNodeWithId( nid_t const parent_id ) {
   forest_[PositionOfNodeInZeroTopology( parent_id )].Coarse( parent_id );
   coarsenings_since_load_balance_++;
}

/**
 * @brief Determines the rank which holds the node of given id.
 * @param id The unique id of the node.
 * @return Rank id for the requested node.
 * @note This function favors feature-envy implementations. It should not be used and rather be a private function.
 */
int TopologyManager::GetRankOfNode( nid_t const id ) const {
   return forest_[PositionOfNodeInZeroTopology( id )].GetRank( id );
}

/**
 * @brief Gives a list which indicates which node should go from which mpi rank onto which mpi rank.
 * @param number_of_ranks The number of ranks available to distribute the load onto.
 * @return A vector of all nodes and their current as well as their future mpi rank.
 */
std::vector<std::tuple<nid_t const, int const, int const>> TopologyManager::PrepareLoadBalancedTopology( int const number_of_ranks ) {

   AssignTargetRankToLeaves( number_of_ranks );
   AssignBalancedLoad();
   std::vector<std::tuple<nid_t const, int const, int const>> ids_current_future_rank_map;

   ListNodeToBalance( ids_current_future_rank_map );

   return ids_current_future_rank_map;
}

/**
 * @brief Indicates whether a node exists in the global Tree, does not make implications about local tree
 * @param id id of the node one is looking for
 * @return true if node exists, false otherwise
 * @note This function favors feature-envy implementations. It should not be used and rather be a private function.
 */
bool TopologyManager::NodeExists( nid_t const id ) const {

   int position_in_zero_topology = PositionOfNodeInZeroTopology( id );
   // Forest size cannot exceed 128x128x128 = 10^6, fits in int ( 10^9 positive values ), cast is safe.
   if( position_in_zero_topology >= 0 && position_in_zero_topology < ( static_cast<int>( forest_.size() ) ) ) {
      return forest_[position_in_zero_topology].NodeExists( id );
   } else {
      return false;
   }
}

/**
 * @brief Gives the current maximum level of any global node.
 *        Can be less than the user set maximum level ( if no interesting physics are present, or at initialization ).
 * @return Globally Maximal Present Level.
 */
unsigned int TopologyManager::GetCurrentMaximumLevel() const {

   std::vector<unsigned int> tree_depths;
   tree_depths.reserve( forest_.size() );
   for( TopologyNode const& tree : forest_ ) {
      tree_depths.emplace_back( tree.GetDepth() );
   }

   return *std::max_element( tree_depths.begin(), tree_depths.end() ) - 1;//Tree Depth starts with 1 by definition, vs levels starts at 0.
}

/**
 * @brief Gives the ids of all globally existing nodes which descent from the given id.
 *        I.e. children, grand-children great-grand-children, ...
 * @param id Unique id of node whose descendants are searched for.
 * @return All ids of globally existing descendants.
 */
std::vector<nid_t> TopologyManager::DescendantIdsOfNode( nid_t const id ) const {

   std::vector<nid_t> descendants;
   std::vector<nid_t> append_list;

   for( auto const& child_id : IdsOfChildren( id ) ) {
      if( NodeExists( child_id ) ) {
         append_list = DescendantIdsOfNode( child_id );
         descendants.insert( descendants.end(), append_list.begin(), append_list.end() );
         descendants.push_back( child_id );
      }
   }

   return descendants;
}

/**
 * @brief Indicates whether the invoking MPI rank holds the given node. $Throws exception if node does not exist!$
 * @param id Unique id of the node to be checked for self-ownership.
 * @param rank The rank on which existence is to be checked.
 * @return true if node is on the same rank as invoker, false otherwise.
 * @note This function favors feature-envy implementations. It should not be used and rather be a private function.
 */
bool TopologyManager::NodeIsOnRank( nid_t const id, int const rank ) const {

#ifndef PERFORMANCE
   if( !NodeExists( id ) ) {
      throw std::logic_error( "Node Ownership cannot be checked - Node does not exist" );
   }
#endif

   return GetRankOfNode( id ) == rank;
}

/**
 * @brief Indicates whether the given node is a leaf or not. $Throws exception if node does not exist!$
 * @param id Unique id of the node to be checked.
 * @return true if node is a leaf, false otherwise.
 */
bool TopologyManager::NodeIsLeaf( nid_t const id ) const {

#ifndef PERFORMANCE
   if( !NodeExists( id ) ) {
      throw std::logic_error( "Node Leaf status cannot be checked - Node does not exist" );
   }
#endif

   return forest_[PositionOfNodeInZeroTopology( id )].NodeIsLeaf( id );
}

/**
 * @brief Determines if the specified node is facing a jump at the given location.
 *        I.e. face does not have a global boundary and no neighbor on the same level exists.
 * @param id Unique id of the node under consideration.
 * @param location Location of interest.
 * @return True if the Face is a Jump, false otherwise.
 */
bool TopologyManager::FaceIsJump( nid_t const id, BoundaryLocation const location ) const {

   if( IsExternalTopologyBoundary( location, id ) ) {
      return false;
   }

   // If the neighbor does not exist and it is not an external BC we have a jump
   return !NodeExists( GetTopologyNeighborId( id, location ) );
}

/**
 * @brief Gives a list of all leaves on this MPI rank
 * @return Local leaf ids.
 */
std::vector<nid_t> TopologyManager::LocalLeafIds() const {
   std::vector<nid_t> local_leaves;
   int const rank_id = MpiUtilities::MyRankId();
   for( TopologyNode const& node : forest_ ) {
      node.LocalLeaves( local_leaves, rank_id );
   }
   return local_leaves;
}

std::vector<nid_t> TopologyManager::LocalInterfaceLeafIds() const {
   std::vector<nid_t> local_leaves;
   int const rank_id = MpiUtilities::MyRankId();
   for( TopologyNode const& node : forest_ ) {
      node.LocalInterfaceLeaves( local_leaves, rank_id );
   }
   return local_leaves;
}

/**
 * @brief Gives a list of all nodes on this MPI rank
 * @return Local node ids.
 */
std::vector<nid_t> TopologyManager::LocalNodeIds() const {
   std::vector<nid_t> local_nodes;
   int const rank_id = MpiUtilities::MyRankId();
   for( TopologyNode const& node : forest_ ) {
      node.LocalNodes( local_nodes, rank_id );
   }
   return local_nodes;
}

/**
 * @brief Gives a list of ids of all globally present leaves.
 * @return Leaf Ids.
 */
std::vector<nid_t> TopologyManager::LeafIds() const {
   std::vector<nid_t> leaves;
   for( TopologyNode const& node : forest_ ) {
      node.GetLeafIds( leaves );
   }
   return leaves;
}

/**
 * @brief Gives the ids of all locally present leaves on the specified level.
 * @param level The level of interest.
 * @return The list of leaf ids.
 */
std::vector<nid_t> TopologyManager::LocalLeafIdsOnLevel( unsigned int const level ) const {
   std::vector<nid_t> leaves;
   int const rank_id = MpiUtilities::MyRankId();
   for( TopologyNode const& node : forest_ ) {
      node.LocalLeavesOnLevel( leaves, rank_id, level );
   }
   return leaves;
}

/**
 * @brief Gives the ids of all globally present leaves on the specified level.
 * @param level The level of interest.
 * @return The list of leaf ids.
 */
std::vector<nid_t> TopologyManager::LeafIdsOnLevel( unsigned int const level ) const {
   std::vector<nid_t> leaves;
   for( TopologyNode const& node : forest_ ) {
      node.GetLeafIdsOnLevel( leaves, level );
   }
   return leaves;
}

void TopologyManager::AssignTargetRanksToLeavesInList( std::vector<nid_t> const& leaves, int const number_of_ranks ) {
   auto const elements_per_rank = ElementsPerRank( leaves.size(), number_of_ranks );
   std::size_t start            = 0;
   for( int rank_id = 0; rank_id < number_of_ranks; ++rank_id ) {
      for( std::size_t i = start; i < start + elements_per_rank[rank_id]; ++i ) {
         forest_[PositionOfNodeInZeroTopology( leaves[i] )].AssignTargetRankToLeaf( leaves[i], rank_id );
      }
      start += elements_per_rank[rank_id];
   }
}

/**
 * @brief Assigns the target rank ( rank on which the node SHOULD reside ) to all leaf nodes.
 *        Uses either a linear or Hilbert Traversal to determine the target rank.
 * @param number_of_ranks The number of ranks available to distribute the load onto.
 */
void TopologyManager::AssignTargetRankToLeaves( int const number_of_ranks ) {

   for( unsigned int level = 0; level <= maximum_level_; ++level ) {
      std::vector<nid_t> leaves = LeafIdsOnLevel( level );
      auto start_multi          = std::partition( std::begin( leaves ), std::end( leaves ), [this]( nid_t const node ) { return forest_[PositionOfNodeInZeroTopology( node )].GetMaterials( node ).size() == 1; } );
      std::vector<nid_t> multiphase_leaves( start_multi, std::end( leaves ) );
      leaves.erase( start_multi, std::end( leaves ) );
      auto start_levelset = std::partition( std::begin( multiphase_leaves ), std::end( multiphase_leaves ), [this]( nid_t const node ) { return LevelOfNode( node ) < maximum_level_; } );
      std::vector<nid_t> levelset_leaves( start_levelset, std::end( multiphase_leaves ) );
      multiphase_leaves.erase( start_levelset, std::end( multiphase_leaves ) );
      OrderNodeIdsBySpaceFillingCurve( leaves );
      OrderNodeIdsBySpaceFillingCurve( multiphase_leaves );
      OrderNodeIdsBySpaceFillingCurve( levelset_leaves );
      AssignTargetRanksToLeavesInList( leaves, number_of_ranks );
      AssignTargetRanksToLeavesInList( multiphase_leaves, number_of_ranks );
      AssignTargetRanksToLeavesInList( levelset_leaves, number_of_ranks );
   }
}

/**
 * @brief Gives the position of the zero level ancestor node identified by the given id in the forest.
 * @param id The id of the node whose ancestor's position is to be determined
 * @return The index of the ancestor node in the zero topology. -1 If no such node could be found.
 */
int TopologyManager::PositionOfNodeInZeroTopology( nid_t const id ) const {

   nid_t level_zero_id = id;
   while( LevelOfNode( level_zero_id ) != 0 ) {
      level_zero_id = ParentIdOfNode( level_zero_id );
   }

   auto node_iterator = std::find_if( forest_.begin(), forest_.end(), [&level_zero_id]( TopologyNode const& node ) { return node.Id() == level_zero_id; } );

   if( node_iterator == forest_.end() ) {
      return -1;
   } else {
      return std::distance( forest_.begin(), node_iterator );
   }
}

/**
 * @brief Calculates a balanced distribution of nodes among the MPI ranks and assigns the determined rank to the nodes.
 *        Does not directly shift nodes among ranks! "Prepares for sending Load".
 */
void TopologyManager::AssignBalancedLoad() {
   for( TopologyNode& node : forest_ ) {
      node.BalanceTargetRanks();
   }
}

/**
 * @brief Gives a list of all nodes, that need to be balanced, i.e. shifted to another MPI rank.
 * @param ids_current_future_rank_map Indirect return parameter.
 * @note Lists the ranks to be balanced as tuple of their id, their current rank and the rank they are supposed to be shifted to
 */
void TopologyManager::ListNodeToBalance( std::vector<std::tuple<nid_t const, int const, int const>>& ids_current_future_rank_map ) {
   for( TopologyNode& node : forest_ ) {
      node.ListUnbalancedNodes( ids_current_future_rank_map );
   }
}

/**
 * @brief Gives some statistics about the distribution of leaves on the ranks.
 * @param number_of_ranks The number of ranks for which distribution is to be documented.
 * @return Formatted string of the leaf-rank distribution.
 */
std::string TopologyManager::LeafRankDistribution( int const number_of_ranks ) {

   std::string leaf_rank_distribution;
   leaf_rank_distribution.reserve( 31 + maximum_level_ * ( 10 + number_of_ranks * 20 ) );

   std::vector<std::vector<int>> leaves_per_level_per_rank( maximum_level_ + 1, std::vector<int>( number_of_ranks, 0 ) );
   for( unsigned int level = 0; level < maximum_level_ + 1; level++ ) {
      std::vector<nid_t> leaves = LeafIdsOnLevel( level );
      for( nid_t id : leaves ) {
         int const rank = GetRankOfNode( id );
         leaves_per_level_per_rank[level][rank]++;
      }
   }

   leaf_rank_distribution.append( "+++ leave rank distribution +++ " );
   for( unsigned int level = 0; level < maximum_level_ + 1; level++ ) {
      leaf_rank_distribution.append( "Level: " + std::to_string( level ) + " " );
      for( int rank = 0; rank < number_of_ranks; rank++ ) {
         leaf_rank_distribution.append( "Rank: " + std::to_string( rank ) + " --> " + std::to_string( leaves_per_level_per_rank[level][rank] ) + " | " );
      }
      leaf_rank_distribution.append( " - " );
   }

   return leaf_rank_distribution;
}

/**
 * @brief Gives out the ids of all globally existent nodes on the specified level.
 * @param level Level of interest.
 * @return Ids of Nodes on level.
 */
std::vector<nid_t> TopologyManager::GlobalIdsOnLevel( unsigned int const level ) const {
   std::vector<nid_t> ids;
   for( TopologyNode const& node : forest_ ) {
      node.IdsOnLevel( level, ids );
   }
   return ids;
}

/**
 * @brief Gives out the ids of only locally existent nodes on the specifed level for a given rank
 * @param level Level of interest.
 * @param rank_id The rank for which the node ids should be given
 * @return Ids of local Nodes on level.
 */
std::vector<nid_t> TopologyManager::IdsOnLevelOfRank( unsigned int const level, int const rank_id ) const {
   std::vector<nid_t> ids;
   for( TopologyNode const& node : forest_ ) {
      node.LocalIdsOnLevel( level, ids, rank_id );
   }
   return ids;
}

/**
 * @brief Gives whether the node with the given id is a multi-phase node, i.e. contains more than one material.
 * @param id Id of the node in question.
 * @return True if the node is multi-phase, false if it is single-phase.
 * @note This is function differs from querying presence of a levelset. Here also nodes on levels other than maximum may return true.
 *       This function favors feature-envy implementations. It should not be used and rather be a private function.
 */
bool TopologyManager::IsNodeMultiPhase( nid_t const id ) const {
   return forest_[PositionOfNodeInZeroTopology( id )].GetMaterials( id ).size() > 1;
}

/**
 * @brief Adds the given material to the node with the given id.
 * @param id Id of the node the material should be added to.
 * @param material The material to be added to the node.
 */
void TopologyManager::AddMaterialToNode( nid_t const id, MaterialName const material ) {
   std::get<0>( local_added_materials_list_ ).push_back( id );
   std::get<1>( local_added_materials_list_ ).push_back( material );
}

/**
 * @brief Removes the given material from the node with the given id.
 * @param id Id of the node the material should be removed from.
 * @param material The material to be removed from the node.
 */
void TopologyManager::RemoveMaterialFromNode( nid_t const id, MaterialName const material ) {
   std::get<0>( local_removed_materials_list_ ).push_back( id );
   std::get<1>( local_removed_materials_list_ ).push_back( material );
}

/**
 * @brief Gives the sorted materials list of the phases present in the given node.
 * @param id Id of the node in question.
 * @return Vector of the materials in the node.
 */
std::vector<MaterialName> TopologyManager::GetMaterialsOfNode( nid_t const id ) const {
   return forest_[PositionOfNodeInZeroTopology( id )].GetMaterials( id );
}

/**
 * @brief Gives the material in a single phase node.
 * @param id Node id.
 * @return The material.
 */
MaterialName TopologyManager::SingleMaterialOfNode( nid_t const id ) const {
   return forest_[PositionOfNodeInZeroTopology( id )].GetSingleMaterial( id );
}

/**
 * @brief Gives whether the given material is present in the given node.
 * @param node_id Id of the node in question.
 * @param material Material in question.
 * @return True if the material is present in the node, false otherwise.
 */
bool TopologyManager::NodeContainsMaterial( nid_t const node_id, MaterialName const material ) const {
   auto materials      = GetMaterialsOfNode( node_id );
   auto block_iterator = std::find( materials.begin(), materials.end(), material );
   return block_iterator == materials.end() ? false : true;
}

/**
 * @brief Indicates - based on the number of mesh changes since last load balancing - whether or not load balancing should be executed.
 * @return Indicator for load balancing.
 */
bool TopologyManager::IsLoadBalancingNecessary() {
   // Check whether load balancing is required based on CC chosen value
   if( coarsenings_since_load_balance_ >= CC::TCULB() || refinements_since_load_balance_ >= CC::TCULB() ) {
      coarsenings_since_load_balance_ = 0;
      refinements_since_load_balance_ = 0;
      return true;
   } else {
      return false;
   }
}

/**
 * @brief Gives the number of global nodes and leaves in a std::pair
 * @return std::pair<#Nodes, #Leaves>
 */
std::pair<unsigned int, unsigned int> TopologyManager::NodeAndLeafCount() const {
   std::pair<unsigned int, unsigned int> node_leaf_count = std::make_pair<unsigned int, unsigned int>( 0, 0 );
   for( TopologyNode const& node : forest_ ) {
      node.NodeLeafCount( node_leaf_count );
   }
   return node_leaf_count;
}

/**
 * @brief Gives the number of leaves which contain an interface.
 * @return number of interface containing leaves.
 */
unsigned int TopologyManager::InterfaceLeafCount() const {
   return std::accumulate( std::cbegin( forest_ ),
                           std::cend( forest_ ),
                           0u,
                           []( auto const partial_sum, auto const& node ) { return partial_sum + node.InterfaceLeafCount(); } );
}

/**
 * @brief Gives a list of pairs. An entry at index i corresponds to the MPI rank_id i. It lists the node and leaf count on this rank
 * @return Vector of std::pair<#Nodes, #Leaves> of size total number of ranks.
 */
std::vector<std::pair<unsigned int, unsigned int>> TopologyManager::NodesAndLeavesPerRank() const {
   std::vector<std::pair<unsigned int, unsigned int>> nodes_and_leaves_per_rank;
   for( TopologyNode const& node : forest_ ) {
      node.RankWiseNodeLeafCount( nodes_and_leaves_per_rank );
   }
   return nodes_and_leaves_per_rank;
}

/**
 * @brief Gives the number of leaves which contain an interface for each rank.
 * @retrun Vector of interface leaf counts. Postion refelects rank id.
 */
std::vector<unsigned int> TopologyManager::InterfaceLeavesPerRank() const {
   std::vector<unsigned int> interface_leaves_per_rank;
   for( TopologyNode const& node : forest_ ) {
      node.RankWiseInterfaceLeafCount( interface_leaves_per_rank );
   }
   return interface_leaves_per_rank;
}

/**
 * @brief returns the number of nodes and blocks in a std::pair
 * @return std::pair<#Nodes, #Blocks>
 */
std::pair<unsigned int, unsigned int> TopologyManager::NodeAndBlockCount() const {
   std::pair<unsigned int, unsigned int> node_block_count = std::make_pair<unsigned int, unsigned int>( 0, 0 );
   for( TopologyNode const& node : forest_ ) {
      node.NodeBlockCount( node_block_count );
   }
   return node_block_count;
}

/**
 * @brief Gives a list of pairs. An entry at index i corresponds to the MPI rank_id i. It lists the node and block count on this rank
 * std::pair<#Nodes, #Blocks>.
 * @return Vector of std::pair<#Nodes, #Blocks> of size total number of ranks.
 */
std::vector<std::pair<unsigned int, unsigned int>> TopologyManager::NodesAndBlocksPerRank() const {
   std::vector<std::pair<unsigned int, unsigned int>> nodes_and_blocks_per_rank;
   for( TopologyNode const& node : forest_ ) {
      node.RankWiseNodeBlockCount( nodes_and_blocks_per_rank );
   }
   return nodes_and_blocks_per_rank;
}

/**
 * @brief Gives the global count of nodes holding more than one phase in the topology.
 * @return Number of Multiphase nodes
 */
unsigned int TopologyManager::MultiPhaseNodeCount() const {
   unsigned int count = 0;
   for( TopologyNode const& node : forest_ ) {
      count += node.MultiPhaseNodeCount();
   }
   return count;
}

/**
 * @brief Restores the complete topology based on a list of node ids, the number of phases foreach node and the material identifiers of the
 *  respective phases. The topology is also load balanced.
 * @param ids A global list of all nodes that are part of this topology.
 * @param number_of_phases The number of phases present in each node. Has to have the same length as ids.
 * @param materials The material identifiers for all phases. The length equals the accumulation of all entries in number_of_phases.
 * @return A list identifying the nodes that are handled by the current rank by means of their indices in the input list ids.
 */
std::vector<unsigned int> TopologyManager::RestoreTopology( std::vector<nid_t> ids, std::vector<unsigned short> number_of_phases,
                                                            std::vector<unsigned short> materials ) {
   std::array<std::vector<unsigned int>, CC::AMNL()> indices_on_level;
   for( unsigned int index = 0; index < ids.size(); ++index ) {
      indices_on_level[LevelOfNode( ids[index] )].push_back( index );
   }

   // sanity check on level 0 (no PERFORMANCE macro required here since only used during restart of simulations)
   if( indices_on_level[0].size() != forest_.size() ) {
      throw std::runtime_error( "Level-zero topology of input and restart file do not match! (size)" );
   }
   for( auto const index_node : indices_on_level[0] ) {
      if( PositionOfNodeInZeroTopology( ids[index_node] ) == -1 ) {
         throw std::runtime_error( "Level-zero topology of input and restart file do not match! (topology)" );
      }
   }

   // build up the topology tree
   for( unsigned int level = 0; level < CC::AMNL(); ++level ) {
      for( auto const index_node : indices_on_level[level] ) {
         nid_t const id     = ids[index_node];
         TopologyNode& root = forest_[PositionOfNodeInZeroTopology( id )];
         // check whether the parent has to be refined
         if( !root.NodeExists( id ) ) {
            root.Refine( ParentIdOfNode( id ) );// safe to call on level 0 because the nodes always exist
         }
         // assign the node's materials
         unsigned int offset_material = std::accumulate( number_of_phases.begin(), number_of_phases.begin() + index_node, 0 );
         for( unsigned int index_material = offset_material; index_material < offset_material + number_of_phases[index_node]; ++index_material ) {
            root.AddMaterial( id, static_cast<MaterialName>( materials[index_material] ) );
         }
      }
   }

   // load balance topology
   PrepareLoadBalancedTopology( MpiUtilities::NumberOfRanks() );

   // return the indices in the input list of the nodes that ended up on this rank
   std::vector<unsigned int> local_indices;
   int const rank_id = MpiUtilities::MyRankId();
   for( unsigned int index = 0; index < ids.size(); ++index ) {
      if( GetRankOfNode( ids[index] ) == rank_id ) {
         local_indices.push_back( index );
      }
   }
   return local_indices;
}

/**
 * @brief Returns a list of all neighboring leaf ids of a node in a given direction.
 * @param node_id Id of the node in question.
 * @param direction Direction in which to check for neighbors.
 * @return List of tuple of node ids of neighboring leaves and their tree level difference compared to node_id
 */
std::vector<nid_t> TopologyManager::GetNeighboringLeaves( nid_t const node_id, BoundaryLocation const direction ) const {

   std::vector<nid_t> id_list;
   if( !IsExternalBoundary( direction, node_id, number_of_nodes_on_level_zero_ ) ) {
      nid_t neighbor_id = GetNeighborId( node_id, direction );
      if( NodeExists( neighbor_id ) ) {
         //find set of lowest level neighbors ( leaves )
         std::function<bool( nid_t const )> sibling_function;
         switch( direction ) {
            case BoundaryLocation::West: {
               sibling_function = EastInSiblingPack;
            } break;
            case BoundaryLocation::East: {
               sibling_function = WestInSiblingPack;
            } break;
            case BoundaryLocation::North: {
               sibling_function = SouthInSiblingPack;
            } break;
            case BoundaryLocation::South: {
               sibling_function = NorthInSiblingPack;
            } break;
            case BoundaryLocation::Bottom: {
               sibling_function = TopInSiblingPack;
            } break;
#ifndef PERFORMANCE
            case BoundaryLocation::Top: {
               sibling_function = BottomInSiblingPack;
            } break;
            default:
               throw std::invalid_argument( "Got invalid direction for TopologyManager::GetNeighborLeavesAndLevelDifference()" );
#else
            default: /* BoundaryLocation::Top */ {
               sibling_function = BottomInSiblingPack;
            }
#endif
         }
         std::vector<nid_t> open_neighbor_ids;
         open_neighbor_ids.push_back( neighbor_id );
         while( open_neighbor_ids.size() > 0 ) {
            nid_t open_id = open_neighbor_ids.back();
            open_neighbor_ids.pop_back();
            if( NodeIsLeaf( open_id ) ) {
               //open_id is leaf node -> add to list
               id_list.push_back( open_id );
            } else {
               //open_id has children -> get the relevant ones
               std::vector<nid_t> const open_children_ids = IdsOfChildren( open_id );
               for( nid_t open_children_id : open_children_ids ) {
                  if( sibling_function( open_children_id ) ) {
                     open_neighbor_ids.push_back( open_children_id );
                  }
               }
            }
         }
      } else {
         //neighbor has lower resolution
         while( !NodeExists( neighbor_id ) && neighbor_id > 2 ) {
            neighbor_id = ParentIdOfNode( neighbor_id );
         }
         //should not need test, since node is not supposed to be boundary
         if( neighbor_id > 2 ) id_list.push_back( neighbor_id );
      }
   }
   return id_list;
   //Note: This function could also be implemented recursively over topology nodes
}

/**
 * @brief Gives a rank specific offset, i. e. a count of how many leafs are on lower (by rank id) rank.
 *        (e.g., three ranks with three leaves each. Offset rank 0 = 0, Offset rank 1 = 3, Offset rank 2 = 6)
 * @param rank The rank for which the offset is to be obtained.
 * @return The offset.
 */
unsigned long long int TopologyManager::LeafOffsetOfRank( int const rank ) const {
   std::vector<std::pair<unsigned int, unsigned int>>&& rank_node_map = NodesAndLeavesPerRank();
   auto const cend                                                    = rank_node_map.size() > std::size_t( rank ) ? rank_node_map.cbegin() + rank : rank_node_map.cend();
   return std::accumulate( rank_node_map.cbegin(), cend, 0ll,
                           []( unsigned int const& a, std::pair<unsigned int, unsigned int> const& b ) { return a + b.second; } );
}

/**
 * @brief Gives a rank specific offset, i. e. a count of how many interface leafs are on lower ( by rank id ) rank.
 * @param rank The rank for which the offset is to be obtained.
 * @return The offset.
 */
unsigned long long int TopologyManager::InterfaceLeafOffsetOfRank( int const rank ) const {
   std::vector<unsigned int> rank_node_map = InterfaceLeavesPerRank();
   auto const cend                         = rank_node_map.size() > std::size_t( rank ) ? rank_node_map.cbegin() + rank : rank_node_map.cend();
   return std::accumulate( rank_node_map.cbegin(), cend, 0ll );
}

/**
 * @brief Gives a rank specific offset, i. e. a count of how many nodes are on lower ( by rank id ) rank.
 * @param rank The rank for which the offset is to be obtained.
 * @return The offset.
 */
unsigned long long int TopologyManager::NodeOffsetOfRank( int const rank ) const {
   std::vector<std::pair<unsigned int, unsigned int>>&& rank_node_map = NodesAndLeavesPerRank();
   auto const cend                                                    = rank_node_map.size() > std::size_t( rank ) ? rank_node_map.cbegin() + rank : rank_node_map.cend();
   return std::accumulate( rank_node_map.cbegin(), cend, 0ll,
                           []( unsigned int const& a, std::pair<unsigned int, unsigned int> const& b ) { return a + b.first; } );
}

/**
 * @brief Gives a rank specific offset, i. e. a count of how many nodes and block are on lower ( by rank id ) rank.
 * @param rank The rank for which the offset is to be obtained.
 * @return The offset.
 */
std::pair<unsigned long long int, unsigned long long int> TopologyManager::NodeAndBlockOffsetOfRank( int const rank ) const {
   std::vector<std::pair<unsigned int, unsigned int>>&& rank_node_block_map = NodesAndBlocksPerRank();
   auto const cend                                                          = rank_node_block_map.size() > std::size_t( rank ) ? rank_node_block_map.cbegin() + rank : rank_node_block_map.cend();
   return std::make_pair( std::accumulate( rank_node_block_map.cbegin(), cend, 0u,
                                           []( unsigned int const& a, std::pair<unsigned int, unsigned int> const& b ) { return a + b.first; } ),
                          std::accumulate( rank_node_block_map.cbegin(), cend, 0u,
                                           []( unsigned int const& a, std::pair<unsigned int, unsigned int> const& b ) { return a + b.second; } ) );
}

/**
    * @brief Determines whether a location ( including edges and corners ) of a block is at the edge of the computational domain.
    * @param location The  direction of the edge under consideration.
    * @param id The id of the node under investigation.
    * @return True if the edge is a domain edge, false otherwise, i.e. internal edge.
    * @note Does not check for dimensionality! I. e. callers responsibility to only call on existing locations ( e. g. NOT Top in 1D ).
    */
bool TopologyManager::IsExternalTopologyBoundary( BoundaryLocation const location, nid_t const id ) const {
   return PeriodicIsExternalBoundary( location, id, number_of_nodes_on_level_zero_, active_periodic_locations_ );
}

/**
    * @brief Gives the maximum level
    * @return Maximum level
    */
unsigned int TopologyManager::GetMaximumLevel() const {
   return maximum_level_;
}
/**
    * @brief Gives the number of nodes on level zero
    * @return Number of nodes
    */
std::array<unsigned int, 3> TopologyManager::GetNumberOfNodesOnLevelZero() const {
   return number_of_nodes_on_level_zero_;
}

/**
   * @brief Gives the id of a neighbor at the provided direction.
   * @param id The id of the node whose neighbor is to be found.
   * @param location Direction in which the neighbor is located.
   * @return Id of the neighbor.
   */
nid_t TopologyManager::GetTopologyNeighborId( nid_t const id, BoundaryLocation const location ) const {
   return GetPeriodicNeighborId( id, location, number_of_nodes_on_level_zero_, active_periodic_locations_ );
}