//===----------------------------- tree.h ---------------------------------===//
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
#ifndef TREE_H
#define TREE_H

#include "block_definitions/interface_block.h"
#include "node.h"
#include "topology_manager.h"
#include <memory>
#include <unordered_map>
#include <vector>

/**
 * @brief The Tree class holds the information about local material data. Data
 * on the MR levels is stored in Node containers as Binary-, Quad-, Oct-tree in
 * one, two and three dimensions, respectively. Tree does not manipulate the
 * data itself, it is an accessor that provides data to a solver or similar. No
 * real tree-searches are performed due to the unique node indexing. Tree is
 * rather a guarantee to access the correct data and that changes are
 * transferable e.g. to the TopologyManager in the proper way. The tree must not
 * change global data at any time.
 */
class Tree {
  // topology for global information of the nodes on all ranks
  TopologyManager const &topology_;
  // size of a cubic block on level zero (non-dimensionalized)
  double const node_size_on_level_zero_;
  // all nodes contained in this tree ( current rank )
  std::vector<std::unordered_map<nid_t, Node>> nodes_;

  void InsertNode(nid_t const id, std::vector<MaterialName> const materials,
                  std::int8_t const interface_tag);

public:
  Tree() = delete;
  explicit Tree(TopologyManager const &topology,
                unsigned int const maximum_level,
                double const node_size_on_level_zero);
  ~Tree() = default;
  Tree(Tree const &) = delete;
  Tree &operator=(Tree const &) = delete;
  Tree(Tree &&) = delete;
  Tree &operator=(Tree &&) = delete;

  /**
   * @brief Gives the non-dimensionalized node size on level zero
   * @return Number of nodes
   */
  inline double GetNodeSizeOnLevelZero() const {
    return node_size_on_level_zero_;
  }

  // Functions for manipulating the tree
  Node &CreateNode(
      nid_t const id, std::vector<MaterialName> const &materials,
      std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
      std::unique_ptr<InterfaceBlock> interface_block = nullptr);
  Node &CreateNode(nid_t const id, std::vector<MaterialName> const &materials);
  void RemoveNodeWithId(nid_t const id);
  std::vector<nid_t> RefineNode(nid_t const id);

  // Functions to return leaf nodes
  std::vector<std::reference_wrapper<Node>> Leaves();
  std::vector<std::reference_wrapper<Node const>> Leaves() const;
  std::vector<std::reference_wrapper<Node>>
  LeavesOnLevel(unsigned int const level);
  std::vector<std::reference_wrapper<Node const>>
  LeavesOnLevel(unsigned int const level) const;
  std::vector<std::reference_wrapper<Node>>
  NonLevelsetLeaves(unsigned int const level);
  std::vector<std::reference_wrapper<Node const>>
  NonLevelsetLeaves(unsigned int const level) const;
  std::vector<std::reference_wrapper<Node const>> InterfaceLeaves() const;

  // Functions to return all nodes ( not only leaf nodes )
  std::vector<std::reference_wrapper<Node>>
  NodesOnLevel(unsigned int const level);
  std::vector<std::reference_wrapper<Node const>>
  NodesOnLevel(unsigned int const level) const;
  std::vector<std::reference_wrapper<Node>> NodesWithLevelset();
  std::vector<std::reference_wrapper<Node const>> NodesWithLevelset() const;
  std::vector<std::reference_wrapper<Node const>> AllNodes() const;

  // Functions to return specific ndes
  Node const &GetNodeWithId(nid_t const id) const;
  Node &GetNodeWithId(nid_t const id);
  std::pair<nid_t const, Node> &NodeIdPair(nid_t const id);
  std::pair<nid_t const, Node> const &NodeIdPair(nid_t const id) const;
  std::unordered_map<nid_t, Node> &GetLevelContent(unsigned int const level);
  std::unordered_map<nid_t, Node> const &
  GetLevelContent(unsigned int const level) const;

  /**
   * @brief Gives a reference to the complete node list in this tree instance, i
   * e. the complete tree on current MPI rank.
   * @return List of nodes. An array for each level holding arbitrary number of
   * nodes on each level.
   */
  inline std::vector<std::unordered_map<nid_t, Node>> &FullNodeList() {
    return nodes_;
  }
  /**
   * @brief Const overload.
   */
  inline std::vector<std::unordered_map<nid_t, Node>> const &
  FullNodeList() const {
    return nodes_;
  }
};

#endif // TREE_H
