//===--------------- zero_gradient_boundary_condition.h -------------------===//
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
#ifndef ZERO_GRADIENT_BOUNDARY_CONDITION_H
#define ZERO_GRADIENT_BOUNDARY_CONDITION_H

#include "boundary_constants.h"
#include "levelset_boundary_condition.h"
#include "material_boundary_condition.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief The ZeroGradientBoundaryCondition class implements a zero gradient
 * (extending) external boundary condition of the domain.
 */
template <BoundaryLocation LOC>
class ZeroGradientBoundaryCondition : public MaterialBoundaryCondition,
                                      public LevelsetBoundaryCondition {

  /**
   * @brief Updates the halo cells from the internal cells according to the
   * zero-gradient condition.
   * @param host_buffer Reference of the buffer that is to be updated.
   */
  template <class T>
  inline void
  UpdateZeroGradient(T (&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]) const {
    auto start_indices = BoundaryConstants<LOC>::HaloStartIndices();
    auto end_indices = BoundaryConstants<LOC>::HaloEndIndices();

    for (unsigned int i = start_indices[0]; i < end_indices[0]; ++i) {
      for (unsigned int j = start_indices[1]; j < end_indices[1]; ++j) {
        for (unsigned int k = start_indices[2]; k < end_indices[2]; ++k) {
          host_buffer[i][j][k] =
              BoundaryConstants<LOC>::ZeroGradientValue(host_buffer, i, j, k);
        }
      }
    }
  }

public:
  ZeroGradientBoundaryCondition() = default;
  ~ZeroGradientBoundaryCondition() = default;
  ZeroGradientBoundaryCondition(ZeroGradientBoundaryCondition const &) = delete;
  ZeroGradientBoundaryCondition &
  operator=(ZeroGradientBoundaryCondition const &) = delete;
  ZeroGradientBoundaryCondition(ZeroGradientBoundaryCondition &&) = delete;
  ZeroGradientBoundaryCondition &&
  operator=(ZeroGradientBoundaryCondition &&) = delete;

  /**
   * @brief See base class. Imposes a zero-gradient condition at the boundary.
   */
  void
  UpdateMaterialExternal(Node &node,
                         MaterialFieldType const field_type) const override {
    unsigned int const number_of_fields = MF::ANOF(field_type);
    for (auto &host_mat_block : node.GetPhases()) {
      for (unsigned int field_index = 0; field_index < number_of_fields;
           ++field_index) {
        double(&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] =
            host_mat_block.second.GetFieldBuffer(field_type, field_index);
        UpdateZeroGradient(cells);
      }
    }
  }

  /**
   * @brief See base class. Adjusted to zero-gradient condition.
   */
  void UpdateLevelsetExternal(
      Node &node, InterfaceBlockBufferType const buffer_type) const override {
    /*  NH TODO this if construct should be avoided, therefore different
     * neighbor relations in CommunicationManger needed for levelset vs.
     * Material/Tag Halo updates.
     */
    if (node.HasLevelset()) {
      double(&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetBuffer(buffer_type);
      UpdateZeroGradient(buffer);
    }
  }

  /**
   * @brief See base class. Adjusted to zero-gradient condition.
   */
  void UpdateInterfaceTagExternal(std::int8_t (
      &interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]) const override {
    UpdateZeroGradient(interface_tags);
  }

  /**
   * @brief Identifies the Location of the BoundaryCondition.
   * @return A BoundaryLocation indicating the position of the
   * BoundaryCondition, i. e. which halo cells are updated by it.
   */
  BoundaryLocation GetLocation() const { return LOC; }
};

#endif // ZERO_GRADIENT_BOUNDARY_CONDITION_H
