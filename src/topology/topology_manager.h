//===----------------------- topology_manager.h ---------------------------===//
//
//                                 ALPACA
//
// Part of ALPACA, under the GNU General Public License as published by
// the Free Software Foundation version 3.
// SPDX-License-Identifier: GPL-3.0-only
//
// If using this code in an academic setting, please cite the following:
// @article{hoppe2022parallel,
//  title={A parallel modular computing environment for three-dimensional
//  multiresolution simulations of compressible flows},
//  author={Hoppe, Nils and Adami, Stefan and Adams, Nikolaus A},
//  journal={Computer Methods in Applied Mechanics and Engineering},
//  volume={391},
//  pages={114486},
//  year={2022},
//  publisher={Elsevier}
// }
//
//===----------------------------------------------------------------------===//
#ifndef TOPOLOGY_MANAGER_H
#define TOPOLOGY_MANAGER_H

#include "communication/mpi_utilities.h"
#include "topology/id_periodic_information.h"
#include "topology_node.h"
#include "user_specifications/compile_time_constants.h"
#include <mpi.h>
#include <unordered_map>
#include <vector>

/**
 * @brief The TopologyManager class handles all aspects relevant for MPI (
 * distributed Memory ) parallelization. I.e. overview of data-to-rank maps,
 * sending datatypes, etc. TopologyManager must not trigger changes in local
 * data. It may only be informed about changes in the local trees. In such a
 * case it updates the global information automatically and spreads the
 * information to the threads.
 * @note TopologyManager does not have a link to any tree object as it is not
 * interested in tree details, but only in "node counts".
 */
class TopologyManager {

  unsigned int const maximum_level_;
  unsigned int const active_periodic_locations_;
  std::array<unsigned int, 3> const number_of_nodes_on_level_zero_;

  std::vector<nid_t> local_refine_list_;

  // use tuples of vectors ( instead of the more intuitive vector of tuples ) to
  // ease the use in MPI communication
  std::tuple<std::vector<nid_t>, std::vector<MaterialName>>
      local_added_materials_list_; // List holds id and added materials ( of
                                   // this id ) $USED ONLY IN MULTIPHASE
                                   // VERSION$
  std::tuple<std::vector<nid_t>, std::vector<MaterialName>>
      local_removed_materials_list_; // List holds id and removed materials ( of
                                     // this id ) $USED ONLY IN MULTIPHASE
                                     // VERSION$

  std::unordered_map<nid_t, TopologyNode> forest_;

  unsigned int coarsenings_since_load_balance_;
  unsigned int refinements_since_load_balance_;

  void SetCurrentRanksAccordingToTargetRanks();
  std::vector<std::tuple<nid_t const, int const, int const>> NodesToBalance();

  void AssignTargetRanksToLeavesInList(std::vector<nid_t> const &leaves,
                                       int const number_of_ranks);
  void AssignTargetRankToLeaves(int const number_of_ranks);

  void AssignTargetRankToParents();

public:
  explicit TopologyManager(std::array<unsigned int, 3> level_zero_blocks = {1,
                                                                            1,
                                                                            1},
                           unsigned int const maximum_level = 0,
                           unsigned int active_periodic_locations = 0);
  ~TopologyManager() = default;
  TopologyManager(TopologyManager const &) = delete;
  TopologyManager &operator=(TopologyManager const &) = delete;
  TopologyManager(TopologyManager &&) = delete;
  TopologyManager &operator=(TopologyManager &&) = delete;

  // Counters:
  unsigned int MultiPhaseNodeCount() const;
  unsigned int InterfaceLeafCount() const;
  std::pair<unsigned int, unsigned int> NodeAndLeafCount() const;
  std::pair<unsigned int, unsigned int> NodeAndBlockCount() const;
  // Rank-wise counters:
  std::vector<std::pair<unsigned int, unsigned int>> NodesAndLeavesPerRank(
      int const number_of_ranks = MpiUtilities::NumberOfRanks()) const;
  std::vector<std::pair<unsigned int, unsigned int>> NodesAndBlocksPerRank(
      int const number_of_ranks = MpiUtilities::NumberOfRanks()) const;
  std::vector<unsigned int> InterfaceLeavesPerRank(
      int const number_of_ranks = MpiUtilities::NumberOfRanks()) const;
  // Offset-counters:
  unsigned long long int LeafOffsetOfRank(
      int const rank,
      int const number_of_ranks = MpiUtilities::NumberOfRanks()) const;
  unsigned long long int InterfaceLeafOffsetOfRank(
      int const rank,
      int const number_of_ranks = MpiUtilities::NumberOfRanks()) const;
  unsigned long long int NodeOffsetOfRank(
      int const rank,
      int const number_of_ranks = MpiUtilities::NumberOfRanks()) const;
  std::pair<unsigned long long int, unsigned long long int>
  NodeAndBlockOffsetOfRank(
      int const rank,
      int const number_of_ranks = MpiUtilities::NumberOfRanks()) const;

  // Single node testers:
  bool NodeExists(nid_t const id) const;
  bool NodeIsOnRank(nid_t const id, int const rank) const;
  bool NodeIsLeaf(nid_t const id) const;
  bool IsNodeMultiPhase(nid_t const id) const;
  int GetRankOfNode(nid_t const id) const;

  bool NodeContainsMaterial(nid_t const node_id,
                            MaterialName const material) const;
  std::vector<MaterialName> GetMaterialsOfNode(nid_t const id) const;
  MaterialName SingleMaterialOfNode(nid_t const id) const;

  // Topological questions:
  bool FaceIsJump(nid_t const id, BoundaryLocation const location) const;
  bool IsExternalTopologyBoundary(BoundaryLocation const location,
                                  nid_t const id) const;
  std::vector<nid_t>
  GetNeighboringLeaves(nid_t const id, BoundaryLocation const location) const;
  nid_t GetTopologyNeighborId(nid_t const id,
                              BoundaryLocation const location) const;

  // Simple Getters:
  unsigned int GetMaximumLevel() const;
  std::array<unsigned int, 3> GetNumberOfNodesOnLevelZero() const;
  unsigned int GetCurrentMaximumLevel() const;
  bool IsLoadBalancingNecessary();

  // Node listings:
  std::vector<nid_t> LocalLeafIds() const;
  std::vector<nid_t> LocalInterfaceLeafIds() const;
  std::vector<nid_t> LeafIds() const;
  std::vector<nid_t> LocalLeafIdsOnLevel(unsigned int const level) const;
  std::vector<nid_t> LeafIdsOnLevel(unsigned int const level) const;
  std::vector<nid_t> DescendantIdsOfNode(nid_t const id) const;
  std::vector<nid_t> LocalIds() const;
  std::vector<nid_t> IdsOnLevel(unsigned int const level) const;
  std::vector<nid_t> LocalIdsOnLevel(unsigned int const level) const;

  // Node and/or topology altering:
  void RefineNodeWithId(nid_t const id);
  void CoarseNodeWithId(nid_t const parent_id);
  void AddMaterialToNode(nid_t const id, MaterialName const material);
  void RemoveMaterialFromNode(nid_t const id, MaterialName const material);

  bool UpdateTopology();
  std::vector<std::tuple<nid_t const, int const, int const>>
  PrepareLoadBalancedTopology(int const number_of_ranks);
  std::vector<unsigned int>
  RestoreTopology(std::vector<nid_t> ids,
                  std::vector<unsigned short> number_of_phases,
                  std::vector<MaterialName> materials);

  // Formatted information
  std::string LeafRankDistribution(int const number_of_ranks);
};

#endif // TOPOLOGY_MANAGER_H
