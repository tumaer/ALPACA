//===--------------------- external_halo_manager.h ------------------------===//
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
#ifndef EXTERNAL_HALO_MANAGER_H
#define EXTERNAL_HALO_MANAGER_H

#include "boundary_condition/boundary_specifications.h"
#include "boundary_condition/levelset_boundary_condition.h"
#include "boundary_condition/material_boundary_condition.h"
#include "topology/node.h"
#include <array>
#include <memory>

/**
 * @brief Container of the external boundaries conditions for materials and
 * interfaces. Furthermore, provides functionality to update halo cells of
 * external boundaries of a single node, based on the provided boundary
 * conditions.
 */
class ExternalHaloManager {

private:
  // arrays with full initialized boundary condition on each location (east,
  // west, north, south, top bottom). If not present it is a nullptr.
  std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6> const
      material_boundary_conditions_;
  std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> const
      levelset_boundary_conditions_;

public:
  ExternalHaloManager() = delete;
  explicit ExternalHaloManager(
      std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6>
          material_boundary_conditions,
      std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6>
          levelset_boundary_conditions);
  ~ExternalHaloManager() = default;
  ExternalHaloManager(ExternalHaloManager const &) = delete;
  ExternalHaloManager &operator=(ExternalHaloManager const &) = delete;
  ExternalHaloManager(ExternalHaloManager &&) = delete;
  ExternalHaloManager &operator=(ExternalHaloManager &&) = delete;

  void UpdateLevelsetExternal(Node &node,
                              InterfaceBlockBufferType const buffer_type,
                              BoundaryLocation const loc) const;
  void UpdateInterfaceTagExternal(
      std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
      BoundaryLocation const loc) const;
  void UpdateMaterialExternal(Node &node, MaterialFieldType const field_type,
                              BoundaryLocation const loc) const;
};

#endif /* EXTERNAL_HALO_MANAGER_H */
