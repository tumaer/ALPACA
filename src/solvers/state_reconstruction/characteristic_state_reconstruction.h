//===--------------- characteristic_state_reconstruction.h ----------------===//
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
#ifndef CHARACTERISTIC_STATE_RECONSTRUCTION_H
#define CHARACTERISTIC_STATE_RECONSTRUCTION_H

#include "block_definitions/field_material_definitions.h"
#include "prime_states/prime_state_handler.h"
#include "solvers/state_reconstruction/state_reconstruction.h"

#include "stencils/stencil_utilities.h"

/**
 * @brief Discretization of the spatial reconstruction scheme using
 * characteristic states (obtained by a characteristic decomposition using Roe
 * eigenvectors) for reconstruction.
 */
class CharacteristicStateReconstruction
    : public StateReconstruction<CharacteristicStateReconstruction> {

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
      double const (&Roe_eigenvectors_left)[MF::ANOE()][MF::ANOE()],
      double const (&Roe_eigenvectors_right)[MF::ANOE()][MF::ANOE()],
      double const cell_size, unsigned int const i, unsigned int const j,
      unsigned int const k) const {

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

    constexpr auto conservative_equation_summation_sequence_ =
        MakeConservativeEquationSummationSequence(
            std::make_index_sequence<MF::ANOE() - DTI(CC::DIM())>{});
    std::array<double, ReconstructionStencil::StencilSize()> u_characteristic;
    std::array<double, MF::ANOE()> characteristic_average_plus;
    std::array<double, MF::ANOE()> characteristic_average_minus;

    for (unsigned int n = 0; n < MF::ANOE();
         ++n) { // n is index of characteristic field (eigenvalue, eigenvector)
      // Characteristic decomposition
      for (unsigned int m = 0; m < ReconstructionStencil::StencilSize(); ++m) {
        u_characteristic[m] = 0.0;
        for (unsigned int const l :
             conservative_equation_summation_sequence_[DTI(
                 DIR)]) { // l is index of conservative equation, iterated in
                          // symmetry-preserving sequence
          u_characteristic[m] +=
              Roe_eigenvectors_left[n][l] *
              block.GetAverageBuffer(MF::ASOE()[l])
                  [i + x_reconstruction_offset *
                           (m - ReconstructionStencil::DownstreamStencilSize())]
                  [j + y_reconstruction_offset *
                           (m - ReconstructionStencil::DownstreamStencilSize())]
                  [k +
                   z_reconstruction_offset *
                       (m - ReconstructionStencil::DownstreamStencilSize())];
        } // L-Loop
      }   // M-Loop

      characteristic_average_minus[n] =
          SU::Reconstruction<ReconstructionStencil, SP::UpwindLeft>(
              u_characteristic, cell_size);
      characteristic_average_plus[n] =
          SU::Reconstruction<ReconstructionStencil, SP::UpwindRight>(
              u_characteristic, cell_size);
    } // N-Loop

    auto const reconstructed_conservatives_minus = TransformToPhysicalSpace(
        characteristic_average_minus, Roe_eigenvectors_right);
    auto const reconstructed_conservatives_plus = TransformToPhysicalSpace(
        characteristic_average_plus, Roe_eigenvectors_right);

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
  CharacteristicStateReconstruction() : StateReconstruction() {}
  ~CharacteristicStateReconstruction() = default;
  CharacteristicStateReconstruction(CharacteristicStateReconstruction const &) =
      delete;
  CharacteristicStateReconstruction &
  operator=(CharacteristicStateReconstruction const &) = delete;
  CharacteristicStateReconstruction(CharacteristicStateReconstruction &&) =
      delete;
  CharacteristicStateReconstruction &
  operator=(CharacteristicStateReconstruction &&) = delete;
};

#endif // CHARACTERISTIC_STATE_RECONSTRUCTION_H
