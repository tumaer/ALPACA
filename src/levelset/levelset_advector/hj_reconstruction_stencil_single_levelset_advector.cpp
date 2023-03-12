//===----- hj_reconstruction_stencil_single_levelset_advector.cpp ---------===//
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
#include "hj_reconstruction_stencil_single_levelset_advector.h"

#include "enums/interface_tag_definition.h"
#include "stencils/stencil_utilities.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief Calculates the right-hand side to solve the level-set advection
 * equation. The advection equation corresponds to equation 12 in \cite Hu2006
 * @note The implemented version differs from \cite Fedkiw2000 since the stencil
 * average is performed only once on the gradient instead of twice on the
 * levelset values.
 * @param node See base class.
 */
void HjReconstructionStencilSingleLevelsetAdvector::AdvectImplementation(
    Node &node) const {

  using ReconstructionStencil = ReconstructionStencilSetup::Concretize<
      levelset_reconstruction_stencil>::type;

  InterfaceBlock &interface_block = node.GetInterfaceBlock();
  double const cell_size = node.GetCellSize();
  double const one_cell_size = 1.0 / cell_size;

  std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
  double const(&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetReinitializedBuffer(InterfaceDescription::Levelset);
  double(&levelset_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetRightHandSideBuffer(InterfaceDescription::Levelset);

  double const(&interface_velocity)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetInterfaceStateBuffer(InterfaceState::Velocity);

  double derivatives[DTI(CC::DIM())][2];
  for (unsigned int d = 0; d < DTI(CC::DIM()); ++d) {
    for (unsigned int e = 0; e < 2; ++e) {
      derivatives[d][e] = 0.0;
    }
  }

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

          double const u_interface =
              interface_velocity[i][j][k] * one_cell_size;

          derivatives[0][0] = SU::Derivative<ReconstructionStencil,
                                             SP::UpwindLeft, Direction::X>(
              levelset_reinitialized, i, j, k, 1.0);
          derivatives[0][1] = SU::Derivative<ReconstructionStencil,
                                             SP::UpwindRight, Direction::X>(
              levelset_reinitialized, i, j, k, 1.0);

          if constexpr (CC::DIM() != Dimension::One) {
            derivatives[1][0] = SU::Derivative<ReconstructionStencil,
                                               SP::UpwindLeft, Direction::Y>(
                levelset_reinitialized, i, j, k, 1.0);
            derivatives[1][1] = SU::Derivative<ReconstructionStencil,
                                               SP::UpwindRight, Direction::Y>(
                levelset_reinitialized, i, j, k, 1.0);
          }

          if constexpr (CC::DIM() == Dimension::Three) {
            derivatives[2][0] = SU::Derivative<ReconstructionStencil,
                                               SP::UpwindLeft, Direction::Z>(
                levelset_reinitialized, i, j, k, 1.0);
            derivatives[2][1] = SU::Derivative<ReconstructionStencil,
                                               SP::UpwindRight, Direction::Z>(
                levelset_reinitialized, i, j, k, 1.0);
          }

          double const old_levelset_sign =
              Signum(levelset_reinitialized[i][j][k]);
          double const godunov_hamiltonian =
              GodunovHamiltonian(derivatives, old_levelset_sign);

          levelset_rhs[i][j][k] = -u_interface * godunov_hamiltonian;
        }
      } // i
    }   // j
  }     // k
}
