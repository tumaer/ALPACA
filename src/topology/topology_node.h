//===----------------------- topology_node.h ------------------------------===//
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
#ifndef TOPOLOGY_NODE_H
#define TOPOLOGY_NODE_H

#include "materials/material_definitions.h"
#include "topology/node_id_type.h"
#include <tuple>
#include <vector>

namespace TopologyNodeConstants {
constexpr int unassigned_rank = -1;
} // namespace TopologyNodeConstants

/**
 * @brief The TopologyNode class organizes the light weight global ( over MPI
 * ranks ) node information in a tree structure. Allowing the TopologyManager
 * efficient searches.
 */
class TopologyNode {

  int current_rank_;
  int target_rank_;
  bool is_leaf_;
  std::vector<MaterialName> materials_;

public:
  explicit TopologyNode(
      int const rank = TopologyNodeConstants::unassigned_rank);
  explicit TopologyNode(
      std::vector<MaterialName> const &material,
      int const rank = TopologyNodeConstants::unassigned_rank);
  ~TopologyNode() = default;
  TopologyNode(TopologyNode const &) = delete;
  TopologyNode &operator=(TopologyNode const &) = delete;
  TopologyNode(TopologyNode &&) = delete;
  TopologyNode &operator=(TopologyNode &&) = delete;

  void AddMaterial(MaterialName const material);
  void RemoveMaterial(MaterialName const material);

  std::vector<MaterialName> Materials() const;
  MaterialName SingleMaterial() const;
  std::size_t NumberOfMaterials() const;

  void MakeParent();
  void MakeLeaf();

  bool IsLeaf() const;

  int Rank() const;
  bool IsOnRank(int const rank) const;
  int TargetRank() const;
  void AssignTargetRank(int const rank);
  void SetCurrentRankAccordingToTargetRank();

  bool IsBalanced() const;
};

#endif // TOPOLOGY_NODE_H
