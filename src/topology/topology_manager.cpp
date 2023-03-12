//===----------------------- topology_manager.cpp -------------------------===//
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
#include "topology_manager.h"

#include <algorithm>
#include <functional>
#include <mpi.h>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "communication/mpi_utilities.h"
#include "topology/id_information.h"
#include "topology/node_id_type.h"
#include "topology/space_filling_curve_order.h"
#include "utilities/container_operations.h"
#include "utilities/string_operations.h"

namespace {
/**
 * @brief Gives a count of elements that should be on each respective rank to
 * obtained a balanced load.
 * @param number_of_elements The amount of elements to distribute accross the
 * ranks.
 * @param number_of_ranks The amount of ranks to distribute the elements onto.
 * @return Vector of size number_of_ranks. Each entry gives the amount of
 * elements the respective rank should hold in a well-balanced scenario.
 */
std::vector<std::size_t> ElementsPerRank(std::size_t const number_of_elements,
                                         int const number_of_ranks) {
  // Cast is save. Negative ranks indicate something is broken already anyways.
  std::size_t rank_count = static_cast<std::size_t>(number_of_ranks);
  std::size_t remainder = number_of_elements % rank_count;
  std::vector<std::size_t> elements_per_rank(rank_count);
  for (std::size_t i = 0; i < rank_count; ++i) {
    std::size_t fraction = 0;
    fraction += number_of_elements / rank_count;
    if (i < remainder) {
      fraction++;
    }
    elements_per_rank[i] = fraction;
  }
  return elements_per_rank;
}

/**
 * @brief Checks if the given node is a multiphase node.
 * @param node Topology node that is to be checked for the multiphase condition.
 * @return True if node is multiphase and false otherwise.
 */
bool IsMultiPhase(TopologyNode const &node) {
  return node.NumberOfMaterials() > 1;
}

/**
 * @brief Creates a vector with elements from max_value in descending order
 * until zero included.
 * @param max_value The highest value present in the resulting vector.
 * @return The vector with descending elements.
 */
std::vector<unsigned int> ElementsDescendingFrom(unsigned int const max_value) {
  std::vector<unsigned int> v(max_value + 1);
  std::generate(std::begin(v), std::end(v),
                [value = max_value]() mutable { return value--; });
  return v;
}

} // namespace

/**
 * @brief Default constructor. Creates the local and global Id list on level
 * zero. Default arguments allow testability.
 * @param maximum_level The possible maximum level in this simulation.
 * @param level_zero_nodes_x, level_zero_nodes_y, level_zero_nodes_z Number of
 * blocks on level zero in the x/y/z-axis extension.
 * @param active_periodic_locations Side of the domain on which periodic
 * boundaries are activated.
 */
TopologyManager::TopologyManager(
    std::array<unsigned int, 3> const level_zero_blocks,
    unsigned int const maximum_level,
    unsigned int const active_periodic_locations)
    : maximum_level_(maximum_level),
      active_periodic_locations_(active_periodic_locations),
      number_of_nodes_on_level_zero_(level_zero_blocks), forest_{},
      coarsenings_since_load_balance_{0}, refinements_since_load_balance_{0} {
  nid_t id = IdSeed();

  std::vector<nid_t> initialization_list;

  /*
   * The Nodes are created in a spatial fashion traversing through X, Y and
   * finally Z. The ids in this traversal are not continuous due to the implicit
   * shadow levels in the tree Therefore some non-straightforward index magic
   * needs to be applied Initialization happens only on Level 0
   */
  for (unsigned int i = 0; i < number_of_nodes_on_level_zero_[2]; ++i) {
    for (unsigned int j = 0; j < number_of_nodes_on_level_zero_[1]; ++j) {
      for (unsigned int k = 0; k < number_of_nodes_on_level_zero_[0]; ++k) {
        forest_.emplace(id, 0);
        initialization_list.push_back(id);
        id = EastNeighborOfNodeWithId(
            initialization_list
                .back()); // find eastern neighbor of the just created Node
      }
      // Index magic to create the correct node once the ( outer ) loop counters
      // are resetted
      id = NorthNeighborOfNodeWithId(
          initialization_list[number_of_nodes_on_level_zero_[0] *
                              (j + number_of_nodes_on_level_zero_[1] * i)]);
    }
    id = TopNeighborOfNodeWithId(
        initialization_list[(number_of_nodes_on_level_zero_[1] *
                             number_of_nodes_on_level_zero_[0] * i)]);
  }

  // Assign correct ranks to nodes
  PrepareLoadBalancedTopology(MpiUtilities::NumberOfRanks());
}

/**
 * @brief Updates the topology based on the recorded refinements, coarsening,
 * weights and material changes.
 * @return True if Communication_managers cache needs to be invalidated.
 * @note The coarsening list is specially guarded, it needs to be flushed before
 * it is considered here.
 */
bool TopologyManager::UpdateTopology() {

  // Flag that specifies whether the communication manager should update its
  // neighbor relations
  bool invalidate_communication_manager_cache = false;
  int const number_of_ranks = MpiUtilities::NumberOfRanks();

  // Tree update
  // refine
  std::vector<nid_t> global_refine_list;
  MpiUtilities::LocalToGlobalData(local_refine_list_, MPI_LONG_LONG_INT,
                                  number_of_ranks, global_refine_list);

  for (auto const &refine_id : global_refine_list) {
    TopologyNode &parent = forest_.at(refine_id);
    parent.MakeParent();
    for (auto child_id : IdsOfChildren(refine_id)) {
      forest_.emplace(std::piecewise_construct, std::make_tuple(child_id),
                      std::make_tuple(parent.Rank()));
    }
  }
  // Invalididate cache if any node has been refined
  if (global_refine_list.size() > 0) {
    invalidate_communication_manager_cache = true;
  }

  local_refine_list_.clear();

  refinements_since_load_balance_ += global_refine_list.size();

  // UPDATE MATERIALS OF NODES
  std::tuple<std::vector<nid_t>, std::vector<MaterialName>>
      global_materials_list;
  MpiUtilities::LocalToGlobalData(std::get<0>(local_added_materials_list_),
                                  MPI_LONG_LONG_INT, number_of_ranks,
                                  std::get<0>(global_materials_list));
  MpiUtilities::LocalToGlobalData(std::get<1>(local_added_materials_list_),
                                  MPI_UNSIGNED_SHORT, number_of_ranks,
                                  std::get<1>(global_materials_list));

#ifndef PERFORMANCE
  if (std::get<0>(global_materials_list).size() !=
      std::get<1>(global_materials_list).size()) {
    throw std::logic_error("Unequally sized material-add-lists encountered");
  }
#endif

  for (unsigned int i = 0; i < std::get<0>(global_materials_list).size(); ++i) {
    forest_.at(std::get<0>(global_materials_list)[i])
        .AddMaterial(std::get<1>(global_materials_list)[i]);
  }

  std::get<0>(local_added_materials_list_).clear();
  std::get<1>(local_added_materials_list_).clear();

  // We reuse the global list
  std::get<0>(global_materials_list).clear();
  std::get<1>(global_materials_list).clear();

  MpiUtilities::LocalToGlobalData(std::get<0>(local_removed_materials_list_),
                                  MPI_LONG_LONG_INT, number_of_ranks,
                                  std::get<0>(global_materials_list));
  MpiUtilities::LocalToGlobalData(std::get<1>(local_removed_materials_list_),
                                  MPI_UNSIGNED_SHORT, number_of_ranks,
                                  std::get<1>(global_materials_list));

#ifndef PERFORMANCE
  if (std::get<0>(global_materials_list).size() !=
      std::get<1>(global_materials_list).size()) {
    throw std::logic_error("Created list of unequal length");
  }
#endif

  for (unsigned int i = 0; i < std::get<0>(global_materials_list).size(); ++i) {
    forest_.at(std::get<0>(global_materials_list)[i])
        .RemoveMaterial(std::get<1>(global_materials_list)[i]);
  }

  std::get<0>(local_removed_materials_list_).clear();
  std::get<1>(local_removed_materials_list_).clear();

  return invalidate_communication_manager_cache;
}

/**
 * @brief Marks the node with the given id for refinement.
 * @param id The id of the leaf that is to be refined.
 * @note The actual refinement happens in a bundled fashion in another function.
 * $No checks for leaves are performed caller is responsible$
 */
void TopologyManager::RefineNodeWithId(nid_t const id) {
  local_refine_list_.push_back(id);
}

/**
 * @brief Adds the parent whose children may be coarsened to the coarsening
 * list. The actual data deletion is than bundled using this list.
 * @param parent_id The id of the node that is to be made a leaf.
 */
void TopologyManager::CoarseNodeWithId(nid_t const parent_id) {
  auto const child_ids(IdsOfChildren(parent_id));
  std::for_each(
      std::cbegin(child_ids), std::cend(child_ids),
      [&forest = forest_](auto const child_id) { forest.erase(child_id); });
  forest_.at(parent_id).MakeLeaf();
  coarsenings_since_load_balance_++;
}

/**
 * @brief Determines the rank which holds the node of given id.
 * @param id The unique id of the node.
 * @return Rank id for the requested node.
 * @note This function favors feature-envy implementations. It should not be
 * used and rather be a private function.
 */
int TopologyManager::GetRankOfNode(nid_t const id) const {
  return forest_.at(id).Rank();
}

/**
 * @brief Gives a list which indicates which node should go from which mpi rank
 * onto which mpi rank.
 * @param number_of_ranks The number of ranks available to distribute the load
 * onto.
 * @return A vector of all nodes and their current as well as their target mpi
 * rank.
 */
std::vector<std::tuple<nid_t const, int const, int const>>
TopologyManager::PrepareLoadBalancedTopology(int const number_of_ranks) {
  AssignTargetRankToLeaves(number_of_ranks);
  AssignTargetRankToParents();
  auto nodes_to_balance = NodesToBalance();
  SetCurrentRanksAccordingToTargetRanks();
  return nodes_to_balance;
}

/**
 * @brief Indicates whether a node exists in the global Tree, does not make
 * implications about local tree
 * @param id id of the node one is looking for
 * @return true if node exists, false otherwise
 * @note This function favors feature-envy implementations. It should not be
 * used and rather be a private function.
 */
bool TopologyManager::NodeExists(nid_t const id) const {
  return forest_.contains(id);
}

/**
 * @brief Gives the current maximum level of any global node.
 *        Can be less than the user set maximum level ( if no interesting
 * physics are present, or at initialization ).
 * @return Globally Maximal Present Level.
 */
unsigned int TopologyManager::GetCurrentMaximumLevel() const {
  auto const it =
      std::max_element(std::cbegin(forest_), std::cend(forest_),
                       [](auto const &id_node_a, auto const &id_node_b) {
                         return std::get<0>(id_node_a) < std::get<0>(id_node_b);
                       });
  return LevelOfNode(std::get<0>(*it));
}

/**
 * @brief Gives the ids of all globally existing nodes which descent from the
 * given id. I.e. children, grand-children great-grand-children, ...
 * @param id Unique id of node whose descendants are searched for.
 * @return All ids of globally existing descendants.
 */
std::vector<nid_t> TopologyManager::DescendantIdsOfNode(nid_t const id) const {

  std::vector<nid_t> descendants;
  std::vector<nid_t> append_list;

  for (auto const &child_id : IdsOfChildren(id)) {
    if (NodeExists(child_id)) {
      append_list = DescendantIdsOfNode(child_id);
      descendants.insert(descendants.end(), append_list.begin(),
                         append_list.end());
      descendants.push_back(child_id);
    }
  }

  return descendants;
}

/**
 * @brief Indicates whether the invoking MPI rank holds the given node. $Throws
 * exception if node does not exist!$
 * @param id Unique id of the node to be checked for self-ownership.
 * @param rank The rank on which existence is to be checked.
 * @return true if node is on the same rank as invoker, false otherwise.
 * @note This function favors feature-envy implementations. It should not be
 * used and rather be a private function.
 */
bool TopologyManager::NodeIsOnRank(nid_t const id, int const rank) const {

#ifndef PERFORMANCE
  if (!NodeExists(id)) {
    throw std::logic_error(
        "Node Ownership cannot be checked - Node does not exist");
  }
#endif

  return forest_.at(id).IsOnRank(rank);
}

/**
 * @brief Indicates whether the given node is a leaf or not. $Throws exception
 * if node does not exist!$
 * @param id Unique id of the node to be checked.
 * @return true if node is a leaf, false otherwise.
 */
bool TopologyManager::NodeIsLeaf(nid_t const id) const {
  return forest_.at(id).IsLeaf();
}

/**
 * @brief Determines if the specified node is facing a jump at the given
 * location. I.e. face does not have a global boundary and no neighbor on the
 * same level exists.
 * @param id Unique id of the node under consideration.
 * @param location Location of interest.
 * @return True if the Face is a Jump, false otherwise.
 */
bool TopologyManager::FaceIsJump(nid_t const id,
                                 BoundaryLocation const location) const {

  if (IsExternalTopologyBoundary(location, id)) {
    return false;
  }

  // If the neighbor does not exist and it is not an external BC we have a jump
  return !NodeExists(GetTopologyNeighborId(id, location));
}

/**
 * @brief Gives a list of all leaves on this MPI rank
 * @return Local leaf ids.
 */
std::vector<nid_t> TopologyManager::LocalLeafIds() const {
  std::vector<nid_t> local_leaves;
  local_leaves.reserve((forest_.size() / MpiUtilities::NumberOfRanks()) +
                       1); //+1 acts as integer-ceil.
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_),
      std::back_inserter(local_leaves),
      [rank_id = MpiUtilities::MyRankId()](auto const &in) {
        return std::get<1>(in).IsLeaf() && std::get<1>(in).IsOnRank(rank_id);
      },
      [](auto const &in) { return std::get<0>(in); });
  return local_leaves;
}

std::vector<nid_t> TopologyManager::LocalInterfaceLeafIds() const {
  std::vector<nid_t> local_interface_leaves;
  local_interface_leaves.reserve(
      (forest_.size() / MpiUtilities::NumberOfRanks()) +
      1); //+1 acts as integer-ceil.( LocalLeafIds() );
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_),
      std::back_inserter(local_interface_leaves),
      [rank_id = MpiUtilities::MyRankId()](auto const &in) {
        auto const &node = std::get<1>(in);
        return IsMultiPhase(node) && node.IsLeaf() && node.IsOnRank(rank_id);
      },
      [](auto const &in) { return std::get<0>(in); });
  return local_interface_leaves;
}

/**
 * @brief Gives a list of all nodes on this MPI rank
 * @return Local node ids.
 */
std::vector<nid_t> TopologyManager::LocalIds() const {
  std::vector<nid_t> local_nodes;
  local_nodes.reserve((forest_.size() / MpiUtilities::NumberOfRanks()) +
                      1); //+1 acts as integer-ceil.
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_), std::back_inserter(local_nodes),
      [rank = MpiUtilities::MyRankId()](auto const &in) {
        return std::get<1>(in).IsOnRank(rank);
      },
      [](auto const &in) { return std::get<0>(in); });
  return local_nodes;
}

/**
 * @brief Gives a list of ids of all globally present leaves.
 * @return Leaf Ids.
 */
std::vector<nid_t> TopologyManager::LeafIds() const {
  std::vector<nid_t> leaves;
  leaves.reserve(forest_.size());
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_), std::back_inserter(leaves),
      [](auto const &in) { return std::get<1>(in).IsLeaf(); },
      [](auto const &in) { return std::get<0>(in); });
  return leaves;
}

/**
 * @brief Gives the ids of all locally present leaves on the specified level.
 * @param level The level of interest.
 * @return The list of leaf ids.
 */
std::vector<nid_t>
TopologyManager::LocalLeafIdsOnLevel(unsigned int const level) const {
  std::vector<nid_t> local_leaves;
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_),
      std::back_inserter(local_leaves),
      [level, rank = MpiUtilities::MyRankId()](auto const &id_node) {
        auto const &[id, node] = id_node;
        return LevelOfNode(id) == level && node.IsLeaf() && node.IsOnRank(rank);
      },
      [](auto const &in) { return std::get<0>(in); });
  return local_leaves;
}

/**
 * @brief Gives the ids of all globally present leaves on the specified level.
 * @param level The level of interest.
 * @return The list of leaf ids.
 */
std::vector<nid_t>
TopologyManager::LeafIdsOnLevel(unsigned int const level) const {
  std::vector<nid_t> leaves;
  leaves.reserve((forest_.size() / MpiUtilities::NumberOfRanks()) +
                 1); //+1 acts as integer-ceil.
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_), std::back_inserter(leaves),
      [level](auto const &id_node) {
        auto const &[id, node] = id_node;
        return LevelOfNode(id) == level && node.IsLeaf();
      },
      [](auto const &in) { return std::get<0>(in); });
  return leaves;
}

/**
 * @brief Gives out the ids of all globally existent nodes on the specified
 * level.
 * @param level Level of interest.
 * @return Ids of Nodes on level.
 */
std::vector<nid_t> TopologyManager::IdsOnLevel(unsigned int const level) const {
  std::vector<nid_t> ids;
  ids.reserve(forest_.size() /
              ((maximum_level_ - level) +
               2)); // The idea is that most ids are on the finest level
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_), std::back_inserter(ids),
      [level](auto const &in) { return LevelOfNode(std::get<0>(in)) == level; },
      [](auto const &in) { return std::get<0>(in); });
  return ids;
}

/**
 * @brief Gives the ids of locally existent nodes on the specifed level for a
 * given rank.
 * @param level Level of interest.
 * @return Ids of local Nodes on level.
 */
std::vector<nid_t>
TopologyManager::LocalIdsOnLevel(unsigned int const level) const {
  std::vector<nid_t> ids;
  ids.reserve((forest_.size() / MpiUtilities::NumberOfRanks()) +
              1); //+1 acts as integer-ceil.
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_), std::back_inserter(ids),
      [level, rank = MpiUtilities::MyRankId()](auto const &id_node) {
        auto const &[id, node] = id_node;
        return LevelOfNode(id) == level && node.IsOnRank(rank);
      },
      [](auto const &in) { return std::get<0>(in); });
  return ids;
}

/**
 * @brief Assigns the target rank to leaves ( rank on which the leaf SHOULD
 * reside ) such that leaves are distributed among all ranks equally.
 * @param leaves The list of leaves that are to be assigned with a target rank.
 * @param number_of_ranks The number of ranks available to distribute the load
 * onto.
 */
void TopologyManager::AssignTargetRanksToLeavesInList(
    std::vector<nid_t> const &leaves, int const number_of_ranks) {
  auto const elements_per_rank =
      ElementsPerRank(leaves.size(), number_of_ranks);
  std::size_t start = 0;
  for (int rank_id = 0; rank_id < number_of_ranks; ++rank_id) {
    for (std::size_t i = start; i < start + elements_per_rank[rank_id]; ++i) {
      forest_.at(leaves[i]).AssignTargetRank(rank_id);
    }
    start += elements_per_rank[rank_id];
  }
}

/**
 * @brief Assigns the target rank ( rank on which the node SHOULD reside ) based
 * on a space-filling curve to all leaf nodes.
 * @param number_of_ranks The number of ranks available to distribute the load
 * onto.
 */
void TopologyManager::AssignTargetRankToLeaves(int const number_of_ranks) {

  for (unsigned int level = 0; level <= maximum_level_; ++level) {
    std::vector<nid_t> leaves = LeafIdsOnLevel(level);
    // On maximum levels all multies are levelset nodes on coarser levels no
    // levelset exists
    auto start_multi =
        std::partition(std::begin(leaves), std::end(leaves),
                       [&forest = forest_](nid_t const node_id) {
                         return !IsMultiPhase(forest.at(node_id));
                       });
    std::vector<nid_t> multiphase_leaves(start_multi, std::end(leaves));
    leaves.erase(start_multi, std::end(leaves));
    OrderNodeIdsBySpaceFillingCurve(leaves);
    OrderNodeIdsBySpaceFillingCurve(multiphase_leaves);
    AssignTargetRanksToLeavesInList(leaves, number_of_ranks);
    AssignTargetRanksToLeavesInList(multiphase_leaves, number_of_ranks);
  }
}

/**
 * @brief Takes the most frequent rank among children nodes and assigns it as
 * the target rank of their parent. This is done on parents of all levels.
 */
void TopologyManager::AssignTargetRankToParents() {
  std::vector<unsigned int> const descending_levels =
      ElementsDescendingFrom(maximum_level_);
  for (auto const level : descending_levels) {
    for (auto &[id, node] : forest_) {
      if (!node.IsLeaf() && LevelOfNode(id) == level) {
        std::unordered_map<int, std::size_t> child_rank_counter;
        child_rank_counter.reserve(CC::NOC());
        for (auto const cid : IdsOfChildren(id)) {
          child_rank_counter[forest_.at(cid).TargetRank()]++;
        }
        node.AssignTargetRank(std::get<0>(*std::max_element(
            std::cbegin(child_rank_counter), std::cend(child_rank_counter),
            [](auto const &a, auto const &b) {
              return std::get<1>(a) < std::get<1>(b);
            })));
      }
    }
  }
}
/**
 * @brief Iterates through the topology and sets the current to match the target
 * rank.
 */
void TopologyManager::SetCurrentRanksAccordingToTargetRanks() {
  std::for_each(std::begin(forest_), std::end(forest_), [](auto &in) {
    std::get<1>(in).SetCurrentRankAccordingToTargetRank();
  });
}

/**
 * @brief Gives a list of all nodes, that need to be balanced, i.e. shifted to
 * another MPI rank.
 * @return List of all nodes to be transferred to another rank.
 * @note Lists the ranks to be balanced as tuple of their id, their current rank
 * and the rank they are supposed to be shifted to
 */
std::vector<std::tuple<nid_t const, int const, int const>>
TopologyManager::NodesToBalance() {
  std::vector<std::tuple<nid_t const, int const, int const>>
      ids_current_target_rank_map;
  ids_current_target_rank_map.reserve(
      (forest_.size() / MpiUtilities::NumberOfRanks()) +
      1); //+1 acts as integer-ceil.
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_),
      std::back_inserter(ids_current_target_rank_map),
      [](auto const &in) { return !std::get<1>(in).IsBalanced(); },
      [](auto const &in) {
        auto const &[id, node] = in;
        return std::make_tuple(id, node.Rank(), node.TargetRank());
      });
  return ids_current_target_rank_map;
}

/**
 * @brief Gives some statistics about the distribution of leaves on the ranks.
 * @param number_of_ranks The number of ranks for which distribution is to be
 * documented.
 * @return Formatted string of the leaf-rank distribution.
 */
std::string TopologyManager::LeafRankDistribution(int const number_of_ranks) {

  std::string leaf_rank_distribution;
  leaf_rank_distribution.reserve(31 +
                                 maximum_level_ * (10 + number_of_ranks * 20));

  std::vector<std::vector<int>> leaves_per_level_per_rank(
      maximum_level_ + 1, std::vector<int>(number_of_ranks, 0));
  for (unsigned int level = 0; level < maximum_level_ + 1; level++) {
    std::vector<nid_t> leaves = LeafIdsOnLevel(level);
    for (nid_t id : leaves) {
      int const rank = GetRankOfNode(id);
      leaves_per_level_per_rank[level][rank]++;
    }
  }

  leaf_rank_distribution.append("+++ leave rank distribution +++ ");
  for (unsigned int level = 0; level < maximum_level_ + 1; level++) {
    leaf_rank_distribution.append("Level: " + std::to_string(level) + " ");
    for (int rank = 0; rank < number_of_ranks; rank++) {
      leaf_rank_distribution.append(
          "Rank: " + std::to_string(rank) + " --> " +
          std::to_string(leaves_per_level_per_rank[level][rank]) + " | ");
    }
    leaf_rank_distribution.append(" - ");
  }

  return leaf_rank_distribution;
}

/**
 * @brief Gives whether the node with the given id is a multi-phase node, i.e.
 * contains more than one material.
 * @param id Id of the node in question.
 * @return True if the node is multi-phase, false if it is single-phase.
 * @note This is function differs from querying presence of a levelset. Here
 * also nodes on levels other than maximum may return true. This function favors
 * feature-envy implementations. It should not be used and rather be a private
 * function.
 */
bool TopologyManager::IsNodeMultiPhase(nid_t const id) const {
  return IsMultiPhase(forest_.at(id));
}

/**
 * @brief Adds the given material to the node with the given id.
 * @param id Id of the node the material should be added to.
 * @param material The material to be added to the node.
 */
void TopologyManager::AddMaterialToNode(nid_t const id,
                                        MaterialName const material) {
  std::get<0>(local_added_materials_list_).push_back(id);
  std::get<1>(local_added_materials_list_).push_back(material);
}

/**
 * @brief Removes the given material from the node with the given id.
 * @param id Id of the node the material should be removed from.
 * @param material The material to be removed from the node.
 */
void TopologyManager::RemoveMaterialFromNode(nid_t const id,
                                             MaterialName const material) {
  std::get<0>(local_removed_materials_list_).push_back(id);
  std::get<1>(local_removed_materials_list_).push_back(material);
}

/**
 * @brief Gives the sorted materials list of the phases present in the given
 * node.
 * @param id Id of the node in question.
 * @return Vector of the materials in the node.
 */
std::vector<MaterialName>
TopologyManager::GetMaterialsOfNode(nid_t const id) const {
  return forest_.at(id).Materials();
}

/**
 * @brief Gives the material in a single phase node.
 * @param id Node id.
 * @return The material.
 */
MaterialName TopologyManager::SingleMaterialOfNode(nid_t const id) const {
  return forest_.at(id).SingleMaterial();
}

/**
 * @brief Gives whether the given material is present in the given node.
 * @param node_id Id of the node in question.
 * @param material Material in question.
 * @return True if the material is present in the node, false otherwise.
 */
bool TopologyManager::NodeContainsMaterial(nid_t const node_id,
                                           MaterialName const material) const {
  auto materials = GetMaterialsOfNode(node_id);
  auto block_iterator = std::find(materials.begin(), materials.end(), material);
  return block_iterator == materials.end() ? false : true;
}

/**
 * @brief Indicates - based on the number of mesh changes since last load
 * balancing - whether or not load balancing should be executed.
 * @return Indicator for load balancing.
 */
bool TopologyManager::IsLoadBalancingNecessary() {
  // Check whether load balancing is required based on CC chosen value
  if (coarsenings_since_load_balance_ >= CC::TCULB() ||
      refinements_since_load_balance_ >= CC::TCULB()) {
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
std::pair<unsigned int, unsigned int>
TopologyManager::NodeAndLeafCount() const {
  unsigned int const node_count = forest_.size();
  unsigned int const leaf_count =
      std::count_if(std::cbegin(forest_), std::cend(forest_),
                    [](auto const &in) { return std::get<1>(in).IsLeaf(); });
  return {node_count, leaf_count};
}

/**
 * @brief Gives the number of leaves which contain an interface.
 * @return number of interface containing leaves.
 */
unsigned int TopologyManager::InterfaceLeafCount() const {
  return std::count_if(
      std::cbegin(forest_), std::cend(forest_), [](auto const &in) {
        return std::get<1>(in).IsLeaf() && IsMultiPhase(std::get<1>(in));
      });
}

/**
 * @brief Gives a list of pairs. An entry at index i corresponds to the MPI
 * rank_id i. It lists the node and leaf count on this rank
 * @param number_of_ranks The number of ranks present in the tree.
 * @return Vector of std::pair<#Nodes, #Leaves> of size total number of ranks.
 */
std::vector<std::pair<unsigned int, unsigned int>>
TopologyManager::NodesAndLeavesPerRank(int const number_of_ranks) const {
  std::vector<std::pair<unsigned int, unsigned int>> nodes_and_leaves_per_rank(
      number_of_ranks, {0, 0});
  for (auto const &[id, node] : forest_) {
    auto &[node_count, leaf_count] = nodes_and_leaves_per_rank.at(node.Rank());
    node_count++;
    if (node.IsLeaf()) {
      leaf_count++;
    }
  }
  return nodes_and_leaves_per_rank;
}

/**
 * @brief Gives the number of leaves which contain an interface for each rank.
 * @param number_of_ranks The number of ranks present in the tree.
 * @return Vector of interface leaf counts. Postion refelects rank id.
 */
std::vector<unsigned int>
TopologyManager::InterfaceLeavesPerRank(int const number_of_ranks) const {
  std::vector<unsigned int> interface_leaves_per_rank(number_of_ranks, 0);
  for (auto const &[id, node] : forest_) {
    if (node.IsLeaf() && IsMultiPhase(node)) {
      interface_leaves_per_rank.at(node.Rank()) += 1;
    }
  }
  return interface_leaves_per_rank;
}

/**
 * @brief returns the number of nodes and blocks in a std::pair
 * @return std::pair<#Nodes, #Blocks>
 */
std::pair<unsigned int, unsigned int>
TopologyManager::NodeAndBlockCount() const {
  unsigned int const node_count = forest_.size();
  unsigned int const block_count =
      std::accumulate(std::cbegin(forest_), std::cend(forest_), 0,
                      [](auto sum, auto const &in) {
                        return sum += std::get<1>(in).NumberOfMaterials();
                      });
  return {node_count, block_count};
}

/**
 * @brief Gives a list of pairs. An entry at index i corresponds to the MPI
 * rank_id i. It lists the node and block count on this rank std::pair<#Nodes,
 * #Blocks>.
 * @param number_of_ranks The number of ranks present in the tree.
 * @return Vector of std::pair<#Nodes, #Blocks> of size total number of ranks.
 */
std::vector<std::pair<unsigned int, unsigned int>>
TopologyManager::NodesAndBlocksPerRank(int const number_of_ranks) const {
  std::vector<std::pair<unsigned int, unsigned int>> nodes_and_blocks_per_rank(
      number_of_ranks, {0, 0});
  for (auto const &[id, node] : forest_) {
    auto &[node_count, block_count] = nodes_and_blocks_per_rank.at(node.Rank());
    node_count++;
    block_count += node.NumberOfMaterials();
  }
  return nodes_and_blocks_per_rank;
}

/**
 * @brief Gives the global count of nodes holding more than one phase in the
 * topology.
 * @return Number of Multiphase nodes
 */
unsigned int TopologyManager::MultiPhaseNodeCount() const {
  return std::count_if(
      std::cbegin(forest_), std::cend(forest_),
      [](auto const &in) { return IsMultiPhase(std::get<1>(in)); });
}

/**
 * @brief Restores the complete topology based on a list of node ids, the number
 * of phases foreach node and the material identifiers of the respective phases.
 * The topology is also load balanced.
 * @param ids A global list of all nodes that are part of this topology.
 * @param number_of_phases The number of phases present in each node. Has to
 * have the same length as ids.
 * @param materials The material identifiers for all phases. The length equals
 * the accumulation of all entries in number_of_phases.
 * @return A list identifying the nodes that are handled by the current rank by
 * means of their indices in the input list ids.
 */
std::vector<unsigned int>
TopologyManager::RestoreTopology(std::vector<nid_t> ids,
                                 std::vector<unsigned short> number_of_phases,
                                 std::vector<MaterialName> materials) {

  forest_.clear();
  constexpr int initial_rank =
      0; // We first assign all nodes to rank 0, then we balance. This gives
         // consistency between ranks.
  for (std::size_t i = 0, material_counter = 0; i < ids.size(); ++i) {
    std::vector<MaterialName> const materials_in_node(
        std::cbegin(materials) + material_counter,
        std::cbegin(materials) + material_counter + number_of_phases.at(i));
    material_counter += number_of_phases.at(i);
    forest_.emplace(std::piecewise_construct, std::forward_as_tuple(ids.at(i)),
                    std::forward_as_tuple(materials_in_node, initial_rank));
  }

  for (auto &[id, node] : forest_) {
    nid_t const first_child_id = IdsOfChildren(id).front();
    if (forest_.contains(first_child_id)) {
      node.MakeParent();
    }
  }

  PrepareLoadBalancedTopology(MpiUtilities::NumberOfRanks());

  // return the indices in the input list of the nodes that ended up on this
  // rank
  std::vector<nid_t> local_indices;
  local_indices.reserve((forest_.size() / MpiUtilities::NumberOfRanks()) +
                        1); //+1 acts as integer-ceil.
  int const my_rank = MpiUtilities::MyRankId();
  ContainerOperations::transform_if(
      std::cbegin(forest_), std::cend(forest_),
      std::back_inserter(local_indices),
      [my_rank](auto const &in) { return std::get<1>(in).Rank() == my_rank; },
      [](auto const &in) { return std::get<0>(in); });

  std::vector<unsigned int> indices_of_local_nodes(local_indices.size());
  std::transform(std::cbegin(local_indices), std::cend(local_indices),
                 std::begin(indices_of_local_nodes), [&ids](auto const id) {
                   return std::distance(
                       std::cbegin(ids),
                       std::find(std::cbegin(ids), std::cend(ids), id));
                 });
  return indices_of_local_nodes;
}

/**
 * @brief Returns a list of all neighboring leaf ids of a node in a given
 * direction.
 * @param node_id Id of the node in question.
 * @param direction Direction in which to check for neighbors.
 * @return List of tuple of node ids of neighboring leaves and their tree level
 * difference compared to node_id
 */
std::vector<nid_t>
TopologyManager::GetNeighboringLeaves(nid_t const node_id,
                                      BoundaryLocation const direction) const {

  std::vector<nid_t> id_list;
  if (!IsExternalBoundary(direction, node_id, number_of_nodes_on_level_zero_)) {
    nid_t neighbor_id = GetNeighborId(node_id, direction);
    if (NodeExists(neighbor_id)) {
      // find set of lowest level neighbors ( leaves )
      std::function<bool(nid_t const)> sibling_function;
      switch (direction) {
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
        throw std::invalid_argument(
            "Got invalid direction for "
            "TopologyManager::GetNeighborLeavesAndLevelDifference()");
#else
      default: /* BoundaryLocation::Top */ {
        sibling_function = BottomInSiblingPack;
      }
#endif
      }
      std::vector<nid_t> open_neighbor_ids;
      open_neighbor_ids.push_back(neighbor_id);
      while (open_neighbor_ids.size() > 0) {
        nid_t open_id = open_neighbor_ids.back();
        open_neighbor_ids.pop_back();
        if (NodeIsLeaf(open_id)) {
          // open_id is leaf node -> add to list
          id_list.push_back(open_id);
        } else {
          // open_id has children -> get the relevant ones
          std::vector<nid_t> const open_children_ids = IdsOfChildren(open_id);
          for (nid_t open_children_id : open_children_ids) {
            if (sibling_function(open_children_id)) {
              open_neighbor_ids.push_back(open_children_id);
            }
          }
        }
      }
    } else {
      // neighbor has lower resolution
      while (!NodeExists(neighbor_id) && neighbor_id > 2) {
        neighbor_id = ParentIdOfNode(neighbor_id);
      }
      // should not need test, since node is not supposed to be boundary
      if (neighbor_id > 2)
        id_list.push_back(neighbor_id);
    }
  }
  return id_list;
  // Note: This function could also be implemented recursively over topology
  // nodes
}

/**
 * @brief Gives a rank specific offset, i. e. a count of how many leafs are on
 * lower (by rank id) rank. (e.g., three ranks with three leaves each. Offset
 * rank 0 = 0, Offset rank 1 = 3, Offset rank 2 = 6)
 * @param rank The rank for which the offset is to be obtained.
 * @param number_of_ranks The total number of ranks for this distribution.
 * @return The offset.
 */
unsigned long long int
TopologyManager::LeafOffsetOfRank(int const rank,
                                  int const number_of_ranks) const {
  std::vector<std::pair<unsigned int, unsigned int>> &&rank_node_map =
      NodesAndLeavesPerRank(number_of_ranks);
  auto const final_iterator = rank_node_map.size() > std::size_t(rank)
                                  ? std::cbegin(rank_node_map) + rank
                                  : std::cend(rank_node_map);
  return std::accumulate(std::cbegin(rank_node_map), final_iterator, 0ll,
                         [](unsigned int const &a,
                            std::pair<unsigned int, unsigned int> const &b) {
                           return a + b.second;
                         });
}

/**
 * @brief Gives a rank specific offset, i. e. a count of how many interface
 * leafs are on lower ( by rank id ) rank.
 * @param rank The rank for which the offset is to be obtained.
 * @param number_of_ranks The total number of ranks for this distribution.
 * @return The offset.
 */
unsigned long long int
TopologyManager::InterfaceLeafOffsetOfRank(int const rank,
                                           int const number_of_ranks) const {
  std::vector<unsigned int> rank_node_map =
      InterfaceLeavesPerRank(number_of_ranks);
  auto const final_iterator = rank_node_map.size() > std::size_t(rank)
                                  ? std::cbegin(rank_node_map) + rank
                                  : std::cend(rank_node_map);
  return std::accumulate(std::cbegin(rank_node_map), final_iterator, 0ll);
}

/**
 * @brief Gives a rank specific offset, i. e. a count of how many nodes are on
 * lower ( by rank id ) rank.
 * @param rank The rank for which the offset is to be obtained.
 * @param number_of_ranks The total number of ranks for this distribution.
 * @return The offset.
 */
unsigned long long int
TopologyManager::NodeOffsetOfRank(int const rank,
                                  int const number_of_ranks) const {
  std::vector<std::pair<unsigned int, unsigned int>> const rank_node_map =
      NodesAndLeavesPerRank(number_of_ranks);
  auto const final_iterator = rank_node_map.size() > std::size_t(rank)
                                  ? std::cbegin(rank_node_map) + rank
                                  : std::cend(rank_node_map);
  return std::accumulate(std::cbegin(rank_node_map), final_iterator, 0ll,
                         [](unsigned long long int const sum,
                            std::pair<unsigned int, unsigned int> const &b) {
                           return sum + b.first;
                         });
}

/**
 * @brief Gives a rank specific offset, i. e. a count of how many nodes and
 * block are on lower ( by rank id ) rank.
 * @param rank The rank for which the offset is to be obtained.
 * @param number_of_ranks The total number of ranks for this distribution.
 * @return The offset.
 */
std::pair<unsigned long long int, unsigned long long int>
TopologyManager::NodeAndBlockOffsetOfRank(int const rank,
                                          int const number_of_ranks) const {
  std::vector<std::pair<unsigned int, unsigned int>> const rank_node_block_map =
      NodesAndBlocksPerRank(number_of_ranks);
  auto const final_iterator = rank_node_block_map.size() > std::size_t(rank)
                                  ? rank_node_block_map.cbegin() + rank
                                  : rank_node_block_map.cend();
  return std::make_pair(
      std::accumulate(rank_node_block_map.cbegin(), final_iterator, 0u,
                      [](unsigned int const &a,
                         std::pair<unsigned int, unsigned int> const &b) {
                        return a + b.first;
                      }),
      std::accumulate(rank_node_block_map.cbegin(), final_iterator, 0u,
                      [](unsigned int const &a,
                         std::pair<unsigned int, unsigned int> const &b) {
                        return a + b.second;
                      }));
}

/**
 * @brief Determines whether a location ( including edges and corners ) of a
 * block is at the edge of the computational domain.
 * @param location The  direction of the edge under consideration.
 * @param id The id of the node under investigation.
 * @return True if the edge is a domain edge, false otherwise, i.e. internal
 * edge.
 * @note Does not check for dimensionality! I. e. callers responsibility to only
 * call on existing locations ( e. g. NOT Top in 1D ).
 */
bool TopologyManager::IsExternalTopologyBoundary(
    BoundaryLocation const location, nid_t const id) const {
  return PeriodicIsExternalBoundary(
      location, id, number_of_nodes_on_level_zero_, active_periodic_locations_);
}

/**
 * @brief Gives the maximum level
 * @return Maximum level
 */
unsigned int TopologyManager::GetMaximumLevel() const { return maximum_level_; }
/**
 * @brief Gives the number of nodes on level zero
 * @return Number of nodes
 */
std::array<unsigned int, 3>
TopologyManager::GetNumberOfNodesOnLevelZero() const {
  return number_of_nodes_on_level_zero_;
}

/**
 * @brief Gives the id of a neighbor at the provided direction.
 * @param id The id of the node whose neighbor is to be found.
 * @param location Direction in which the neighbor is located.
 * @return Id of the neighbor.
 */
nid_t TopologyManager::GetTopologyNeighborId(
    nid_t const id, BoundaryLocation const location) const {
  return GetPeriodicNeighborId(id, location, number_of_nodes_on_level_zero_,
                               active_periodic_locations_);
}
