//===---------------------- state_reconstruction.h ------------------------===//
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
#ifndef STATE_RECONSTRUCTION_H
#define STATE_RECONSTRUCTION_H

#include "user_specifications/stencil_setup.h"

/**
 * @brief Helper function to create the index sequence used to enforce symmetry
 * while summing up conservative equation contributions in characteristic
 * decomposition.
 * @tparam RemainingIndices Zero-based index sequence representing the
 * non-momentum equations. Are transformed into real equation indices.
 * @return The created index sequence.
 */
template <std::size_t... NonMomentumIndices>
constexpr std::array<std::array<unsigned int, MF::ANOE()>, DTI(CC::DIM())>
MakeConservativeEquationSummationSequence(
    std::index_sequence<NonMomentumIndices...> const) {
#if DIMENSION == 1
  return {{{ETI(Equation::MomentumX),
            ETI(MF::NthNonMomentumEquation(NonMomentumIndices))...}}};
#elif DIMENSION == 2
  return {{{ETI(Equation::MomentumX), ETI(Equation::MomentumY),
            ETI(MF::NthNonMomentumEquation(NonMomentumIndices))...},
           {ETI(Equation::MomentumX), ETI(Equation::MomentumY),
            ETI(MF::NthNonMomentumEquation(NonMomentumIndices))...}}};
#else
  return {{{ETI(Equation::MomentumY), ETI(Equation::MomentumZ),
            ETI(Equation::MomentumX),
            ETI(MF::NthNonMomentumEquation(NonMomentumIndices))...},
           {ETI(Equation::MomentumZ), ETI(Equation::MomentumX),
            ETI(Equation::MomentumY),
            ETI(MF::NthNonMomentumEquation(NonMomentumIndices))...},
           {ETI(Equation::MomentumX), ETI(Equation::MomentumY),
            ETI(Equation::MomentumZ),
            ETI(MF::NthNonMomentumEquation(NonMomentumIndices))...}}};
#endif
}

/**
 * @brief Interface to apply the spatial reconsturction scheme.
 */
template <typename DerivedStateReconstruction> class StateReconstruction {

  friend DerivedStateReconstruction;

  explicit constexpr StateReconstruction() = default;

public:
  ~StateReconstruction() = default;
  StateReconstruction(StateReconstruction const &) = delete;
  StateReconstruction &operator=(StateReconstruction const &) = delete;
  StateReconstruction(StateReconstruction &&) = delete;
  StateReconstruction &operator=(StateReconstruction &&) = delete;

  /**
   * @brief Solves the first-order Riemann problem using left and right state
   * vectors.
   * @tparam direction.
   * @tparam reconstruction scheme.
   * @param block block.
   * @param eos applied equation of state.
   * @param Roe_eigenvectors_left, Roe_eigenvectors_right Left and right Roe
   * eigenvector matrices for characteristic decomposition
   * @param cell_size cell size.
   * @param i,j,k cell indizes in block.
   * @return The reconstructed left/right conservative and primitive states.
   */
  template <Direction DIR, ReconstructionStencils RECON>
  std::tuple<std::array<double, MF::ANOE()>, std::array<double, MF::ANOE()>,
             std::array<double, MF::ANOP()>, std::array<double, MF::ANOP()>>
  SolveStateReconstruction(
      Block const &block, EquationOfState const &eos,
      double const (&Roe_eigenvectors_left)[MF::ANOE()][MF::ANOE()],
      double const (&Roe_eigenvectors_right)[MF::ANOE()][MF::ANOE()],
      double const cell_size, unsigned int const i, unsigned int const j,
      unsigned int const k) const {
    return static_cast<DerivedStateReconstruction const &>(*this)
        .template SolveStateReconstructionImplementation<DIR, RECON>(
            block, eos, Roe_eigenvectors_left, Roe_eigenvectors_right,
            cell_size, i, j, k);
  }
};

#endif // STATE_RECONSTRUCTION_H
