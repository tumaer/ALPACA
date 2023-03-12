//===------------------ levelset_boundary_condition.h ---------------------===//
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
#ifndef LEVELSET_BOUNDARY_CONDITION_H
#define LEVELSET_BOUNDARY_CONDITION_H

#include "block_definitions/interface_block.h"
#include "topology/node.h"

/**
 * @brief The LevelsetBoundaryCondition class defines an interface for the
 * boundaries of nodes. Sets boundary conditions into the respective node's halo
 * cells. Boundaries work on the external halos of a node. However, the internal
 * cells and the halo cells are stored in one continuous buffer. This means just
 * certain entries (indices) of this continuous buffer are effected by
 * LevelsetBoundaryCondition classes.
 */
class LevelsetBoundaryCondition {

public:
  LevelsetBoundaryCondition() = default;
  virtual ~LevelsetBoundaryCondition() = default;
  LevelsetBoundaryCondition(LevelsetBoundaryCondition const &) = delete;
  LevelsetBoundaryCondition &
  operator=(LevelsetBoundaryCondition const &) = delete;
  LevelsetBoundaryCondition(LevelsetBoundaryCondition &&) = delete;
  LevelsetBoundaryCondition &operator=(LevelsetBoundaryCondition &&) = delete;

  /**
   * @brief Performs all levelset halo updates in external boundaries.
   * @param node Node on which the halo update is done.
   * @param buffer_type Identifier of the buffer type of the levelset block to
   * be updated.
   */
  virtual void
  UpdateLevelsetExternal(Node &node,
                         InterfaceBlockBufferType const buffer_type) const = 0;

  /**
   * @brief Performs all interface tag halo updates in external boundaries.
   * @param interface_tags buffer for which the halo update is done (indirect
   * return).
   */
  virtual void UpdateInterfaceTagExternal(
      std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]) const = 0;
};

#endif // LEVELSET_BOUNDARY_CONDITION_H
