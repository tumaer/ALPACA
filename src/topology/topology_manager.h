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
#ifndef TOPOLOGY_MANAGER_H
#define TOPOLOGY_MANAGER_H

#include <vector>
#include <mpi.h>
#include "topology_node.h"
#include "user_specifications/compile_time_constants.h"
#include "topology/id_periodic_information.h"

/**
 * @brief The TopologyManager class handles all aspects relevant for MPI ( distributed Memory ) parallelization. I.e. overview of data-to-rank maps, sending datatypes, etc.
 *        TopologyManager must not trigger changes in local data. It may only be informed about changes in the local trees. In such a case
 *        it updates the global information automatically and spreads the information to the threads.
 * @note TopologyManager does not have a link to any tree object as it is not interested in tree details, but only in "node counts".
 */
class TopologyManager {

   unsigned int const maximum_level_;
   unsigned int const active_periodic_locations_;
   std::array<unsigned int, 3> const number_of_nodes_on_level_zero_;

   std::vector<nid_t> local_refine_list_;

   // use tuples of vectors ( instead of the more intuitive vector of tuples ) to ease the use in MPI communication
   std::tuple<std::vector<nid_t>, std::vector<MaterialName>> local_added_materials_list_;  //List holds id and added materials ( of this id ) $USED ONLY IN MULTIPHASE VERSION$
   std::tuple<std::vector<nid_t>, std::vector<MaterialName>> local_removed_materials_list_;//List holds id and removed materials ( of this id ) $USED ONLY IN MULTIPHASE VERSION$

   std::vector<TopologyNode> forest_;// A collection of ( root ) trees is a forest

   unsigned int coarsenings_since_load_balance_;
   unsigned int refinements_since_load_balance_;

   int PositionOfNodeInZeroTopology( nid_t const id ) const;
   void AssignBalancedLoad();
   void ListNodeToBalance( std::vector<std::tuple<nid_t const, int const, int const>>& ids_current_future_rank_map );

   void AssignTargetRanksToLeavesInList( std::vector<nid_t> const& leaves, int const number_of_ranks );
   void AssignTargetRankToLeaves( int const number_of_ranks );

public:
   explicit TopologyManager( std::array<unsigned int, 3> level_zero_blocks = { 1, 1, 1 }, unsigned int const maximum_level = 0, unsigned int active_periodic_locations = 0 );
   ~TopologyManager()                        = default;
   TopologyManager( TopologyManager const& ) = delete;
   TopologyManager& operator=( TopologyManager const& ) = delete;
   TopologyManager( TopologyManager&& )                 = delete;
   TopologyManager& operator=( TopologyManager&& ) = delete;

   // Counters:
   unsigned int MultiPhaseNodeCount() const;
   unsigned int InterfaceLeafCount() const;
   std::pair<unsigned int, unsigned int> NodeAndLeafCount() const;
   std::pair<unsigned int, unsigned int> NodeAndBlockCount() const;
   // Rank-wise counters:
   std::vector<std::pair<unsigned int, unsigned int>> NodesAndLeavesPerRank() const;
   std::vector<std::pair<unsigned int, unsigned int>> NodesAndBlocksPerRank() const;
   std::vector<unsigned int> InterfaceLeavesPerRank() const;
   // Offset-counters:
   unsigned long long int LeafOffsetOfRank( int const rank ) const;
   unsigned long long int InterfaceLeafOffsetOfRank( int const rank ) const;
   unsigned long long int NodeOffsetOfRank( int const rank ) const;
   std::pair<unsigned long long int, unsigned long long int> NodeAndBlockOffsetOfRank( int const rank ) const;

   //Single node testers:
   bool NodeExists( nid_t const id ) const;
   bool NodeIsOnRank( nid_t const id, int const rank ) const;
   bool NodeIsLeaf( nid_t const id ) const;
   bool IsNodeMultiPhase( nid_t const id ) const;
   int GetRankOfNode( nid_t const id ) const;

   bool NodeContainsMaterial( nid_t const node_id, MaterialName const material ) const;
   std::vector<MaterialName> GetMaterialsOfNode( nid_t const id ) const;
   MaterialName SingleMaterialOfNode( nid_t const id ) const;

   // Topological questions:
   bool FaceIsJump( nid_t const id, BoundaryLocation const location ) const;
   bool IsExternalTopologyBoundary( BoundaryLocation const location, nid_t const id ) const;
   std::vector<nid_t> GetNeighboringLeaves( nid_t const id, BoundaryLocation const location ) const;
   nid_t GetTopologyNeighborId( nid_t const id, BoundaryLocation const location ) const;

   // Simple Getters:
   unsigned int GetMaximumLevel() const;
   std::array<unsigned int, 3> GetNumberOfNodesOnLevelZero() const;
   unsigned int GetCurrentMaximumLevel() const;
   bool IsLoadBalancingNecessary();

   // Node listings:
   std::vector<nid_t> LocalLeafIds() const;
   std::vector<nid_t> LocalInterfaceLeafIds() const;
   std::vector<nid_t> LeafIds() const;
   std::vector<nid_t> LocalLeafIdsOnLevel( unsigned int const level ) const;
   std::vector<nid_t> LeafIdsOnLevel( unsigned int const level ) const;
   std::vector<nid_t> DescendantIdsOfNode( nid_t const id ) const;
   std::vector<nid_t> LocalNodeIds() const;
   std::vector<nid_t> GlobalIdsOnLevel( unsigned int const level ) const;
   std::vector<nid_t> IdsOnLevelOfRank( unsigned int const level, int const rank_id ) const;

   // Node and/or topology altering:
   void RefineNodeWithId( nid_t const id );
   void CoarseNodeWithId( nid_t const parent_id );
   void AddMaterialToNode( nid_t const id, MaterialName const material );
   void RemoveMaterialFromNode( nid_t const id, MaterialName const material );

   bool UpdateTopology();
   std::vector<std::tuple<nid_t const, int const, int const>> PrepareLoadBalancedTopology( int const number_of_ranks );
   std::vector<unsigned int> RestoreTopology( std::vector<nid_t> ids, std::vector<unsigned short> number_of_phases, std::vector<unsigned short> materials );

   // Formatted information
   std::string LeafRankDistribution( int const number_of_ranks );
};

#endif// TOPOLOGY_MANAGER_H
