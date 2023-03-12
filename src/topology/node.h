//===----------------------------- node.h ---------------------------------===//
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
#ifndef NODE_H
#define NODE_H

#include <memory>
#include <unordered_map>
#include <vector>

#include "block_definitions/block.h"
#include "block_definitions/interface_block.h"
#include "boundary_condition/boundary_specifications.h"
#include "enums/interface_tag_definition.h"
#include "materials/material_definitions.h"
#include "topology/id_information.h"

/**
 * @brief Nodes are the members in the tree. A node holds a block for every
 * phase it contains; the Block then holds the material data. Node is a
 * container that gathers information common for all phases at a given position,
 * as e.g. Boundary Condition types. Every node has a unique index for
 * identification, in particular with respect to MPI.
 */
class Node {

  double const node_size_;
  std::tuple<double const, double const, double const> const node_coordinates_;
  std::unordered_map<MaterialName, Block> phases_;

  // type std::int8_t due to definition of enum InterfaceTag. Needs to be
  // changed in case the enum type changes.
  std::int8_t interface_tags_[CC::TCX()][CC::TCY()][CC::TCZ()];
  std::int8_t integrated_interface_tags_[CC::TCX()][CC::TCY()][CC::TCZ()];

  std::unique_ptr<InterfaceBlock> interface_block_;

public:
  Node() = delete;
  explicit Node(nid_t const id, double const node_size_on_level_zero,
                std::vector<MaterialName> const materials,
                std::int8_t const initial_interface_tag = ITTI(IT::OldCutCell));
  explicit Node(nid_t const id, double const node_size_on_level_zero,
                std::vector<MaterialName> const materials,
                std::int8_t const (
                    &initial_interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
                std::unique_ptr<InterfaceBlock> interface_block = nullptr);
  Node(Node const &) = delete;
  Node &operator=(Node const &) = delete;
  Node(Node &&) = delete;
  Node &operator=(Node &&) = delete;
  ~Node() = default;

  // Functions to return geometry/topology data of the node
  std::tuple<double const, double const, double const>
  GetBlockCoordinates() const;
  double GetBlockSize() const;
  double GetCellSize() const;

  // Functions to get material data of the node
  Block &GetSinglePhase();
  Block const &GetSinglePhase() const;
  Block &GetPhaseByMaterial(MaterialName const material);
  Block const &GetPhaseByMaterial(MaterialName const material) const;
  MaterialName GetSinglePhaseMaterial() const;
  std::vector<MaterialName> GetMaterials() const;
  std::unordered_map<MaterialName, Block> &GetPhases();
  std::unordered_map<MaterialName, Block> const &GetPhases() const;

  void AddPhase(MaterialName const material);
  void RemovePhase(MaterialName const material);
  bool ContainsMaterial(MaterialName const material) const;

  InterfaceBlock &GetInterfaceBlock();
  InterfaceBlock const &GetInterfaceBlock() const;
  void
  SetInterfaceBlock(std::unique_ptr<InterfaceBlock> interface_block = nullptr);
  bool HasLevelset() const;

  std::int8_t GetUniformInterfaceTag() const;
  template <InterfaceDescriptionBufferType C>
  auto GetInterfaceTags() -> std::int8_t (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  template <InterfaceDescriptionBufferType C>
  auto GetInterfaceTags() const -> std::int8_t
      const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

  auto GetInterfaceTags(InterfaceDescriptionBufferType const type)
      -> std::int8_t (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto GetInterfaceTags(InterfaceDescriptionBufferType const type) const
      -> std::int8_t const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
};

#endif // NODE_H
