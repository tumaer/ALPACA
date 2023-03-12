//===----------- weno_iterative_levelset_reinitializer.cpp ----------------===//
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
#include "weno_iterative_levelset_reinitializer.h"
#include "stencils/stencil_utilities.h"
#include "user_specifications/two_phase_constants.h"

/**
 * @brief Computes the advection velocity used for the iterative
 * reinitialization. The advection velocity is a smoothened signum of the
 * levelset. Smoothening is necessary for small absolute level-set values to not
 * reinitialize there.
 * @param levelset The levelset-value.
 * @return The advection velocity.
 */
inline double ComputeAdvectionVelocity(double const levelset) {
  constexpr double epsilon = 1.0;
  return levelset / std::sqrt(levelset * levelset + (epsilon * epsilon));
}

/**
 * @brief The default constructor for a WenoIterativeLevelsetReinitializer.
 * Calls the default constructor of the base class.
 * @param halo_manager See base class.
 */
WenoIterativeLevelsetReinitializer::WenoIterativeLevelsetReinitializer(
    HaloManager &halo_manager)
    : IterativeLevelsetReinitializerBase(halo_manager) {
  // Empty Constructor, besides call of base class constructor.
}

/**
 * @brief Reinitializes the levelset field of a node using a HJ WENO scheme.
 * @param node The node with levelset block which has to be reinitialized.
 * @param levelset_type  Level-set field type that is reinitialized.
 * @param is_last_stage Return whether it's the last RK stage or not.
 * @return The residuum for the current node.
 */
double WenoIterativeLevelsetReinitializer::ReinitializeSingleNodeImplementation(
    Node &node, InterfaceDescriptionBufferType const levelset_type,
    bool const is_last_stage) const {

  using ReconstructionStencil = ReconstructionStencilSetup::Concretize<
      levelset_reconstruction_stencil>::type;

  std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceTags(levelset_type);
  InterfaceBlock &interface_block = node.GetInterfaceBlock();
  double(&levelset_orig)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetInterfaceDescriptionBuffer(
          levelset_type)[InterfaceDescription::Levelset];
  double const(&levelset_0_orig)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetRightHandSideBuffer(InterfaceDescription::Levelset);

  double reinitialization_rhs[CC::TCX()][CC::TCY()][CC::TCZ()];

  double residuum = 0.0;

  double derivatives[DTI(CC::DIM())][2];
  for (unsigned int d = 0; d < DTI(CC::DIM()); ++d) {
    for (unsigned int e = 0; e < 2; ++e) {
      derivatives[d][e] = 0.0;
    }
  }

  for (unsigned int i = 0; i < CC::TCX(); ++i) {
    for (unsigned int j = 0; j < CC::TCY(); ++j) {
      for (unsigned int k = 0; k < CC::TCZ(); ++k) {
        reinitialization_rhs[i][j][k] = 0.0;
      } // k
    }   // j
  }     // i

  for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        if (std::abs(interface_tags[i][j][k]) > ITTI(IT::NewCutCell) ||
            (ReinitializationConstants::ReinitializeCutCells &&
             is_last_stage)) {

          // We normalize the level-set field on the cell size, thus 1.0 is
          // necessary as cell size for stencil evaluation.
          derivatives[0][0] =
              SU::Derivative<ReconstructionStencil, SP::UpwindLeft,
                             Direction::X>(levelset_orig, i, j, k, 1.0);
          derivatives[0][1] =
              SU::Derivative<ReconstructionStencil, SP::UpwindRight,
                             Direction::X>(levelset_orig, i, j, k, 1.0);

          if constexpr (CC::DIM() != Dimension::One) {
            derivatives[1][0] =
                SU::Derivative<ReconstructionStencil, SP::UpwindLeft,
                               Direction::Y>(levelset_orig, i, j, k, 1.0);
            derivatives[1][1] =
                SU::Derivative<ReconstructionStencil, SP::UpwindRight,
                               Direction::Y>(levelset_orig, i, j, k, 1.0);
          }

          if constexpr (CC::DIM() == Dimension::Three) {
            derivatives[2][0] =
                SU::Derivative<ReconstructionStencil, SP::UpwindLeft,
                               Direction::Z>(levelset_orig, i, j, k, 1.0);
            derivatives[2][1] =
                SU::Derivative<ReconstructionStencil, SP::UpwindRight,
                               Direction::Z>(levelset_orig, i, j, k, 1.0);
          }

          double const old_levelset_sign = Signum(levelset_0_orig[i][j][k]);
          double const godunov_hamiltonian =
              GodunovHamiltonian(derivatives, old_levelset_sign);
          double const advection_velocity =
              ComputeAdvectionVelocity(levelset_0_orig[i][j][k]);
          double const increment = ReinitializationConstants::Dtau *
                                   advection_velocity *
                                   (1.0 - godunov_hamiltonian);
          if (ReinitializationConstants::TrackConvergence &&
              std::abs(interface_tags[i][j][k]) <
                  ITTI(IT::ReinitializationBand)) {
            residuum = std::max(residuum, std::abs(increment));
          }

          reinitialization_rhs[i][j][k] = increment;
        }
      } // k
    }   // j
  }     // i

  // add to level-set and also the residuum
  for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        levelset_orig[i][j][k] += reinitialization_rhs[i][j][k];
      } // k
    }   // j
  }     // i

  return residuum;
}
