//===------------------ material_boundary_condition.h ---------------------===//
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
#ifndef MATERIAL_BOUNDARY_CONDITION_H
#define MATERIAL_BOUNDARY_CONDITION_H

#include "boundary_specifications.h"
#include "topology/node.h"

/**
 * @brief The MaterialBoundaryCondition class Defines an Interface for the
 * boundaries of blocks. Sets boundary conditions into the respective block's
 * halo cells. Internal boundaries also inherit from this abstract class.
 * Boundaries work on the Halos of a Block. However, the internal cells and the
 * halo cells are stored in one contigious buffer. This means just certain
 * entries (indices) of this contigious buffer are effected by
 * MaterialBoundaryCondition classes.
 */
class MaterialBoundaryCondition {

public:
  MaterialBoundaryCondition() = default;
  virtual ~MaterialBoundaryCondition() = default;
  MaterialBoundaryCondition(MaterialBoundaryCondition const &) = delete;
  MaterialBoundaryCondition &
  operator=(MaterialBoundaryCondition const &) = delete;
  MaterialBoundaryCondition(MaterialBoundaryCondition &&) = delete;
  MaterialBoundaryCondition &operator=(MaterialBoundaryCondition &&) = delete;

  /**
   * @brief Updates Halo Cells, based on the externalBoundaryType. Updates all
   * Materials.
   * @param Node The node whose halo cells are to be updated.
   * @param field_type Field identifier for material block buffers.
   */
  virtual void
  UpdateMaterialExternal(Node &node,
                         MaterialFieldType const field_type) const = 0;
};

#endif // MATERIAL_BOUNDARY_CONDITION_H
