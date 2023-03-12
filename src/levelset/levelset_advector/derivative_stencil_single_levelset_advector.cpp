//===---------- derivative_stencil_single_levelset_advector.cpp -----------===//
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
#include <numeric>

#include "derivative_stencil_single_levelset_advector.h"
#include "levelset/multi_phase_manager/two_phase_manager.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief Calculates the right-hand side to solve the level-set advection
 * equation with a derivative stencil as proposed in \cite Nourgaliev2007.
 * @param node See base class.
 */
void DerivativeStencilSingleLevelsetAdvector::AdvectImplementation(
    Node &node) const {

  using DerivativeStencil =
      DerivativeStencilSetup::Concretize<derivative_stencil>::type;

  InterfaceBlock &interface_block = node.GetInterfaceBlock();
  double const cell_size = node.GetCellSize();
  double const one_cell_size = 1.0 / cell_size;

  std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
  double const(&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetBaseBuffer(InterfaceDescription::Levelset);
  double const(&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetReinitializedBuffer(InterfaceDescription::Levelset);
  double(&levelset_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetRightHandSideBuffer(InterfaceDescription::Levelset);

  double const(&interface_velocity)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetInterfaceStateBuffer(InterfaceState::Velocity);

  std::array<double, DTI(CC::DIM())> interface_velocity_projection;
  std::array<double, DTI(CC::DIM())> levelset_derivative;
  for (unsigned int d = 0; d < DTI(CC::DIM()); ++d) {
    interface_velocity_projection[d] = 0.0;
    levelset_derivative[d] = 0.0;
  }
  std::array<double, DTI(CC::DIM())> increments;

  for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
        /***
         * In order to calculate the right-hand side for the level-set advection
         * in the 2nd RK stage we need reasonable level-set values in the whole
         * extension band. Thus, it is necessary to calculate level-set
         * advection rhs values for the extension band.
         */
        if (std::abs(interface_tags[i][j][k]) <= ITTI(IT::ExtensionBand)) {

          // Calculate normal to determine x-,y- and z-component of
          // interface_velocity
          std::array<double, 3> const normal =
              GetNormal(levelset_reinitialized, i, j, k);

          double const u_interface =
              interface_velocity[i][j][k] * one_cell_size;

          interface_velocity_projection[0] = u_interface * normal[0];
          levelset_derivative[0] =
              SU::ReconstructionWithUpwinding<DerivativeStencil, Direction::X>(
                  levelset, i, j, k, interface_velocity_projection[0],
                  cell_size);
          increments[0] =
              -interface_velocity_projection[0] * levelset_derivative[0];

          if constexpr (CC::DIM() != Dimension::One) {
            interface_velocity_projection[1] = u_interface * normal[1];
            levelset_derivative[1] =
                SU::ReconstructionWithUpwinding<DerivativeStencil,
                                                Direction::Y>(
                    levelset, i, j, k, interface_velocity_projection[1],
                    cell_size);
            increments[1] =
                -interface_velocity_projection[1] * levelset_derivative[1];
          }

          if constexpr (CC::DIM() == Dimension::Three) {
            interface_velocity_projection[2] = u_interface * normal[2];
            levelset_derivative[2] =
                SU::ReconstructionWithUpwinding<DerivativeStencil,
                                                Direction::Z>(
                    levelset, i, j, k, interface_velocity_projection[2],
                    cell_size);
            increments[2] =
                -interface_velocity_projection[2] * levelset_derivative[2];
          }

          levelset_rhs[i][j][k] = ConsistencyManagedSum(increments);
        }
      } // i
    }   // j
  }     // k
}
