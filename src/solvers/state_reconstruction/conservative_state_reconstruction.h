//===--------------- conservative_state_reconstruction.h ------------------===//
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
#ifndef CONSERVATIVE_STATE_RECONSTRUCTION_H
#define CONSERVATIVE_STATE_RECONSTRUCTION_H

#include "block_definitions/field_material_definitions.h"
#include "prime_states/prime_state_handler.h"
#include "solvers/state_reconstruction/state_reconstruction.h"

#include "stencils/stencil_utilities.h"

/**
 * @brief Discretization of the spatial reconstruction scheme using conservative
 * states for reconstruction.
 */
class ConservativeStateReconstruction
    : public StateReconstruction<ConservativeStateReconstruction> {

  friend StateReconstruction;

  /**
   * @brief Procedure to reconstruct the conservatives/primitive states at cell
   * faces.
   * @tparam DIR spatial direction the reconstruction has to be performed.
   * @param block Block of the phase under consideration.
   * @param eos Underlying equation of state of the phase under consideration
   * used to convert primes and conservatives.
   * @param roe_eigenvectors_left .
   * @param roe_eigenvectors_right .
   * @param cell_size .
   * @param i .
   * @param k .
   * @param j .
   * @return tuple containing left and right reconstructed primitive and
   * conservative states.
   */
  template <Direction DIR, ReconstructionStencils RECON>
  std::tuple<std::array<double, MF::ANOE()>, std::array<double, MF::ANOE()>,
             std::array<double, MF::ANOP()>, std::array<double, MF::ANOP()>>
  SolveStateReconstructionImplementation(
      Block const &block, EquationOfState const &eos,
      double const (&)[MF::ANOE()][MF::ANOE()],
      double const (&)[MF::ANOE()][MF::ANOE()], double const cell_size,
      unsigned int const i, unsigned int const j, unsigned int const k) const {

    using ReconstructionStencil =
        typename ReconstructionStencilSetup::Concretize<RECON>::type;
    constexpr unsigned int x_reconstruction_offset =
        DIR == Direction::X ? 1 : 0;
    constexpr unsigned int y_reconstruction_offset =
        DIR == Direction::Y ? 1 : 0;
    constexpr unsigned int z_reconstruction_offset =
        DIR == Direction::Z ? 1 : 0;
    std::array<double, MF::ANOP()> reconstructed_primes_minus;
    std::array<double, MF::ANOP()> reconstructed_primes_plus;

    std::array<double, ReconstructionStencil::StencilSize()>
        reconstruction_array;
    std::array<double, MF::ANOE()> reconstructed_conservatives_minus;
    std::array<double, MF::ANOE()> reconstructed_conservatives_plus;
    for (unsigned int n = 0; n < MF::ANOE(); ++n) {
      for (unsigned int m = 0; m < ReconstructionStencil::StencilSize(); ++m) {
        reconstruction_array[m] =
            block.GetAverageBuffer(MF::ASOE()[n])
                [i + x_reconstruction_offset *
                         (m - ReconstructionStencil::DownstreamStencilSize())]
                [j + y_reconstruction_offset *
                         (m - ReconstructionStencil::DownstreamStencilSize())]
                [k + z_reconstruction_offset *
                         (m - ReconstructionStencil::DownstreamStencilSize())];
      } // M-Loop
      reconstructed_conservatives_minus[n] =
          SU::Reconstruction<ReconstructionStencil, SP::UpwindLeft>(
              reconstruction_array, cell_size);
      reconstructed_conservatives_plus[n] =
          SU::Reconstruction<ReconstructionStencil, SP::UpwindRight>(
              reconstruction_array, cell_size);
    } // N-Loop
    // To check for invalid cells due to ghost fluid method
    if (reconstructed_conservatives_minus[ETI(Equation::Mass)] <=
            std::numeric_limits<double>::epsilon() ||
        reconstructed_conservatives_plus[ETI(Equation::Mass)] <=
            std::numeric_limits<double>::epsilon())
      return std::make_tuple(
          reconstructed_conservatives_minus, reconstructed_conservatives_plus,
          reconstructed_primes_minus, reconstructed_primes_plus);
    ConservativesToPrimeStates(eos, reconstructed_conservatives_minus,
                               reconstructed_primes_minus);
    ConservativesToPrimeStates(eos, reconstructed_conservatives_plus,
                               reconstructed_primes_plus);
    return std::make_tuple(
        reconstructed_conservatives_minus, reconstructed_conservatives_plus,
        reconstructed_primes_minus, reconstructed_primes_plus);
  }

public:
  ConservativeStateReconstruction() : StateReconstruction() {}
  ~ConservativeStateReconstruction() = default;
  ConservativeStateReconstruction(ConservativeStateReconstruction const &) =
      delete;
  ConservativeStateReconstruction &
  operator=(ConservativeStateReconstruction const &) = delete;
  ConservativeStateReconstruction(ConservativeStateReconstruction &&) = delete;
  ConservativeStateReconstruction &
  operator=(ConservativeStateReconstruction &&) = delete;
};

#endif // CONSERVATIVE_STATE_RECONSTRUCTION_H
