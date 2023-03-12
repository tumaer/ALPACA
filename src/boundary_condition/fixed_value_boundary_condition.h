//===----------------- fixed_value_boundary_condition.h -------------------===//
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
#ifndef FIXED_VALUE_BOUNDARY_CONDITION_H
#define FIXED_VALUE_BOUNDARY_CONDITION_H

#include "boundary_constants.h"
#include "material_boundary_condition.h"

/**
 * @brief The FixedValueBoundaryCondition class imposes a pre-defined value for
 * each variable in all halo cells.
 */
template <BoundaryLocation LOC>
class FixedValueBoundaryCondition : public MaterialBoundaryCondition {

private:
  // fixed values for primestates and conservatives provided at the boundaries
  // (for each material contained in the simulation)
  std::vector<std::array<double, MF::ANOE()>> const fixed_conservatives_;
  std::vector<std::array<double, MF::ANOP()>> const fixed_prime_states_;

public:
  FixedValueBoundaryCondition() = delete;
  /**
   * @brief Default constructor. See base class constructor.
   */
  explicit FixedValueBoundaryCondition(
      std::vector<std::array<double, MF::ANOE()>> const &fixed_conservatives,
      std::vector<std::array<double, MF::ANOP()>> const &fixed_prime_states)
      : fixed_conservatives_(fixed_conservatives),
        fixed_prime_states_(fixed_prime_states) {
    /** Empty besides initializer list */
  }
  ~FixedValueBoundaryCondition() = default;
  FixedValueBoundaryCondition(FixedValueBoundaryCondition const &) = delete;
  FixedValueBoundaryCondition &
  operator=(FixedValueBoundaryCondition const &) = delete;
  FixedValueBoundaryCondition(FixedValueBoundaryCondition &&) = delete;
  FixedValueBoundaryCondition &
  operator=(FixedValueBoundaryCondition &&) = delete;

  /**
   * @brief Imposes predefined values onto the respective halo cells. See base
   * class.
   */
  void
  UpdateMaterialExternal(Node &node,
                         MaterialFieldType const field_type) const override {
    auto start_indices = BoundaryConstants<LOC>::HaloStartIndices();
    auto end_indices = BoundaryConstants<LOC>::HaloEndIndices();

    for (auto &host_mat_block : node.GetPhases()) {
      unsigned int const material_index(MTI(host_mat_block.first));
      switch (field_type) {
      case MaterialFieldType::Conservatives: {
        for (Equation const eq : MF::ASOE()) {
          double(&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] =
              host_mat_block.second.GetRightHandSideBuffer(eq);
          for (unsigned int i = start_indices[0]; i < end_indices[0]; ++i) {
            for (unsigned int j = start_indices[1]; j < end_indices[1]; ++j) {
              for (unsigned int k = start_indices[2]; k < end_indices[2]; ++k) {
                cells[i][j][k] = fixed_conservatives_[material_index][ETI(eq)];
              }
            }
          }
        }
      } break;
      case MaterialFieldType::PrimeStates: {
        for (PrimeState const ps : MF::ASOP()) {
          double(&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] =
              host_mat_block.second.GetPrimeStateBuffer(ps);
          for (unsigned int i = start_indices[0]; i < end_indices[0]; ++i) {
            for (unsigned int j = start_indices[1]; j < end_indices[1]; ++j) {
              for (unsigned int k = start_indices[2]; k < end_indices[2]; ++k) {
                cells[i][j][k] = fixed_prime_states_[material_index][PTI(ps)];
              }
            }
          }
        }
      } break;
      case MaterialFieldType::Parameters: {
        throw std::runtime_error(
            "For the material field parameters fixed value boundary conditions "
            "are not implemented yet!");
      } break;
      default:
        throw std::runtime_error(
            "Material field type not known for fixed value BC!");
      }
    }
  }

  /**
   * @brief Identifies the Location of the BoundaryCondition.
   * @return A BoundaryLocation indicating the position of the
   * BoundaryCondition, i. e. which halo cells are updated by it.
   */
  BoundaryLocation GetLocation() const { return LOC; }
};

#endif // FIXED_VALUE_BOUNDARY_CONDITION_H
