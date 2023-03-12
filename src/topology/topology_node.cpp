//===----------------------- topology_node.cpp ----------------------------===//
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
#include "topology_node.h"

#include "topology/id_information.h"
#include "user_specifications/compile_time_constants.h"
#include "utilities/container_operations.h"
#include <algorithm>

namespace TNC = TopologyNodeConstants;

/**
 * @brief Constructs a topology node as leaf without children, without assigning
 * a future rank and without materials.
 * @param rank Rank holding the node to be created.
 */
TopologyNode::TopologyNode(int const rank)
    : current_rank_(rank), target_rank_(TNC::unassigned_rank), is_leaf_(true),
      materials_() {}

/**
 * @brief Constructs a topology node as leaf with the given materials, but
 * without assigning a future rank.
 * @param materials The materials to be present in the node.
 * @param rank Rank holding the node to be created.
 */
TopologyNode::TopologyNode(std::vector<MaterialName> const &materials,
                           int const rank)
    : current_rank_(rank), target_rank_(TNC::unassigned_rank), is_leaf_(true),
      materials_(ContainerOperations::SortedCopy(materials)) {}

/**
 * @brief Adds the given material to the node.
 * @param material The material to be added to the node.
 * @note Does not perform checks if the material is already present. Inserting
 * the same material twice will lead to undefined behavior.
 */
void TopologyNode::AddMaterial(MaterialName const material) {
  materials_.insert(
      std::upper_bound(std::begin(materials_), std::end(materials_), material),
      material);
}

/**
 * @brief Removes the given material from the node.
 * @param material The material to be removed from the node.
 */
void TopologyNode::RemoveMaterial(MaterialName const material) {
  materials_.erase(
      std::remove(std::begin(materials_), std::end(materials_), material),
      std::end(materials_));
}

/**
 * @brief Gives the  materials present in the node.
 * @return The materials.
 */
std::vector<MaterialName> TopologyNode::Materials() const { return materials_; }

/**
 * @brief Gives the material of a single phase node.
 * @return The material.
 * @note If called on a multi-phase node result is undefined.
 */
MaterialName TopologyNode::SingleMaterial() const { return materials_.front(); }

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
void TopologyNode::MakeParent() { is_leaf_ = false; }

/**
 * @brief Converts a parent node into leaf.
 */
void TopologyNode::MakeLeaf() { is_leaf_ = true; }

/**
 * @brief Whether the node is a leaf.
 * @return True if the node is a leaf, false otherwise.
 */
bool TopologyNode::IsLeaf() const { return is_leaf_; }

/**
 * @brief Gives the rank the node resides on.
 * @return The rank id.
 */
int TopologyNode::Rank() const { return current_rank_; }

/**
 * @brief Indicates if the node resides on the provided rank.
 * @return True if node is on the rank. False otherwise.
 */
bool TopologyNode::IsOnRank(int const rank) const {
  return current_rank_ == rank;
}

/**
 * @brief Gives the rank on which the node should resides on (in the future).
 * @return The rank id.
 */
int TopologyNode::TargetRank() const { return target_rank_; }

/**
 * @brief Sets the desired rank of the node to the given rank.
 * @param rank The rank the node should reside on.
 */
void TopologyNode::AssignTargetRank(int const rank) { target_rank_ = rank; }

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
bool TopologyNode::IsBalanced() const { return current_rank_ == target_rank_; }
