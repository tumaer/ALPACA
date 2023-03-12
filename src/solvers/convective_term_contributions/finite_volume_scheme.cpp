//===--------------------- finite_volume_scheme.cpp -----------------------===//
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
#include "solvers/convective_term_contributions/finite_volume_scheme.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief Standard constructor using an already existing MaterialManager and
 * EigenDecomposition object.
 * @param material_manager .
 * @param eigendecomposition_calculator .
 */
FiniteVolumeScheme::FiniteVolumeScheme(
    MaterialManager const &material_manager,
    EigenDecomposition const &eigendecomposition_calculator)
    : ConvectiveTermSolver(material_manager, eigendecomposition_calculator),
      riemann_solver_(material_manager, eigendecomposition_calculator_),
      state_reconstruction_() {
  /* Empty besides initializer list*/
}

/**
 * @brief Solving the convective term of the system. Using dimension splitting
 * for fluxes in x, y, and z- direction. Also See base class.
 * @note Hotpath function.
 */
void FiniteVolumeScheme::UpdateImplementation(
    std::pair<MaterialName const, Block> const &mat_block,
    double const cell_size,
    double (&fluxes_x)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
    double (&fluxes_y)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
    double (&fluxes_z)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
    double (
        &volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()]) const {

  double roe_eigenvectors_left[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                              [MF::ANOE()][MF::ANOE()];
  double roe_eigenvectors_right[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                               [MF::ANOE()][MF::ANOE()];
  double roe_eigenvalues[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                        [MF::ANOE()];

  double u_hllc_x[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1];
  double u_hllc_y[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1];
  double u_hllc_z[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1];

  constexpr bool require_eigendecomposition =
      (active_equations == EquationSet::NavierStokes ||
       active_equations == EquationSet::Euler) &&
      state_reconstruction_type == StateReconstructionType::Characteristic;

  if constexpr (require_eigendecomposition) {
    // Setting the eigenvectors, eigenvalues to zero is necessary for two-phase
    // simulations as not every entry is necessarily set during
    // eigendecomposition
    for (unsigned int i = 0; i < CC::ICX() + 1; ++i) {
      for (unsigned int j = 0; j < CC::ICY() + 1; ++j) {
        for (unsigned int k = 0; k < CC::ICZ() + 1; ++k) {
          for (unsigned int e = 0; e < MF::ANOE(); ++e) {
            for (unsigned int f = 0; f < MF::ANOE(); ++f) {
              roe_eigenvectors_left[i][j][k][e][f] = 0.0;
              roe_eigenvectors_right[i][j][k][e][f] = 0.0;
            }
            roe_eigenvalues[i][j][k][e] = 0.0;
          }
        }
      }
    }

    eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::X>(
        mat_block, roe_eigenvectors_left, roe_eigenvectors_right,
        roe_eigenvalues);
  }
  ComputeFluxes<Direction::X>(mat_block, fluxes_x, u_hllc_x,
                              roe_eigenvectors_left, roe_eigenvectors_right,
                              cell_size);

  if constexpr (CC::DIM() != Dimension::One) {
    if constexpr (require_eigendecomposition) {
      eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::Y>(
          mat_block, roe_eigenvectors_left, roe_eigenvectors_right,
          roe_eigenvalues);
    }
    ComputeFluxes<Direction::Y>(mat_block, fluxes_y, u_hllc_y,
                                roe_eigenvectors_left, roe_eigenvectors_right,
                                cell_size);
  }

  if constexpr (CC::DIM() == Dimension::Three) {
    if constexpr (require_eigendecomposition) {
      eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::Z>(
          mat_block, roe_eigenvectors_left, roe_eigenvectors_right,
          roe_eigenvalues);
    }
    ComputeFluxes<Direction::Z>(mat_block, fluxes_z, u_hllc_z,
                                roe_eigenvectors_left, roe_eigenvectors_right,
                                cell_size);
  }

  if constexpr (active_equations == EquationSet::GammaModel) {
    // Add source term for non-conservative advections
    for (unsigned int i = 0; i < CC::ICX(); ++i) {
      for (unsigned int j = 0; j < CC::ICY(); ++j) {
        for (unsigned int k = 0; k < CC::ICZ(); ++k) {
          volume_forces[ETI(Equation::Gamma)][i][j][k] -=
              DimensionAwareConsistencyManagedSum(
                  u_hllc_x[i][j + 1][k + 1] - u_hllc_x[i + 1][j + 1][k + 1],
                  u_hllc_y[i + 1][j][k + 1] - u_hllc_y[i + 1][j + 1][k + 1],
                  u_hllc_z[i + 1][j + 1][k] - u_hllc_z[i + 1][j + 1][k + 1]) *
              mat_block.second.GetAverageBuffer(
                  Equation::Gamma)[i + CC::FICX()][j + CC::FICY()]
                                  [k + CC::FICZ()] /
              cell_size;
          volume_forces[ETI(Equation::Pi)][i][j][k] -=
              DimensionAwareConsistencyManagedSum(
                  u_hllc_x[i][j + 1][k + 1] - u_hllc_x[i + 1][j + 1][k + 1],
                  u_hllc_y[i + 1][j][k + 1] - u_hllc_y[i + 1][j + 1][k + 1],
                  u_hllc_z[i + 1][j + 1][k] - u_hllc_z[i + 1][j + 1][k + 1]) *
              mat_block.second.GetAverageBuffer(
                  Equation::Pi)[i + CC::FICX()][j + CC::FICY()]
                               [k + CC::FICZ()] /
              cell_size;
        } // k
      }   // j
    }     // i
  }
}

/**
 * @brief Computes the convective cell face fluxes by performing a finite-volume
 * state reconstruction and solving a Riemann problem afterwards.
 * @param mat_block The block and material information of the phase under
 * consideration.
 * @param fluxes Reference to an array which is filled with the computed fluxes
 * (indirect return parameter).
 * @param roe_eigenvectors_left .
 * @param roe_eigenvectors_right .
 * @param cell_size .
 * @tparam DIR Indicates which spatial direction is to be computed.
 * @note Hotpath function.
 */
template <Direction DIR>
void FiniteVolumeScheme::ComputeFluxes(
    std::pair<MaterialName const, Block> const &mat_block,
    double (&fluxes)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
    double (&u_hllc)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
    double const (
        &Roe_eigenvectors_left)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                               [MF::ANOE()][MF::ANOE()],
    double const (
        &Roe_eigenvectors_right)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                                [MF::ANOE()][MF::ANOE()],
    double const cell_size) const {

  constexpr unsigned int x_start =
      DIR == Direction::X ? CC::FICX() - 1 : CC::FICX();
  constexpr unsigned int y_start =
      DIR == Direction::Y ? CC::FICY() - 1 : CC::FICY();
  constexpr unsigned int z_start =
      DIR == Direction::Z ? CC::FICZ() - 1 : CC::FICZ();

  constexpr unsigned int x_end = CC::LICX();
  constexpr unsigned int y_end = CC::LICY();
  constexpr unsigned int z_end = CC::LICZ();

  constexpr int total_to_internal_offset_x = CC::FICX() - 1;
  constexpr int total_to_internal_offset_y =
      CC::DIM() != Dimension::One ? static_cast<int>(CC::FICY()) - 1 : -1;
  constexpr int total_to_internal_offset_z =
      CC::DIM() == Dimension::Three ? static_cast<int>(CC::FICZ()) - 1 : -1;

  // Access the pair's elements directly.
  auto const &[material, block] = mat_block;

  for (unsigned int i = x_start; i <= x_end; ++i) {
    for (unsigned int j = y_start; j <= y_end; ++j) {
      for (unsigned int k = z_start; k <= z_end; ++k) {
        // Shifted indices to match block index system and roe-ev index system
        int const i_index = i - total_to_internal_offset_x;
        int const j_index = j - total_to_internal_offset_y;
        int const k_index = k - total_to_internal_offset_z;

        auto const [reconstructed_conservatives_left,
                    reconstructed_conservatives_right,
                    reconstructed_primes_left, reconstructed_primes_right] =
            state_reconstruction_.SolveStateReconstruction<
                DIR, reconstruction_stencil>(
                block,
                material_manager_.GetMaterial(material).GetEquationOfState(),
                Roe_eigenvectors_left[i_index][j_index][k_index],
                Roe_eigenvectors_right[i_index][j_index][k_index], cell_size, i,
                j, k);
        // To check for invalid cells due to ghost fluid method
        if constexpr (active_equations != EquationSet::GammaModel) {
          double const B = active_equations == EquationSet::Isentropic
                               ? 0.0
                               : material_manager_.GetMaterial(material)
                                     .GetEquationOfState()
                                     .B();
          if (reconstructed_conservatives_left[ETI(Equation::Mass)] <=
                  std::numeric_limits<double>::epsilon() ||
              reconstructed_conservatives_right[ETI(Equation::Mass)] <=
                  std::numeric_limits<double>::epsilon())
            continue;
          if (reconstructed_primes_left[PTI(PrimeState::Pressure)] <= -B ||
              reconstructed_primes_right[PTI(PrimeState::Pressure)] <= -B)
            continue;
        }

        if constexpr (active_equations == EquationSet::GammaModel) {
          auto const [face_fluxes, u_hllc_single] =
              riemann_solver_.SolveGammaRiemannProblem<DIR>(
                  material, reconstructed_conservatives_left,
                  reconstructed_conservatives_right, reconstructed_primes_left,
                  reconstructed_primes_right);
          for (unsigned int n = 0; n < MF::ANOE(); ++n) {
            fluxes[n][i_index][j_index][k_index] += face_fluxes[n];
          }
          u_hllc[i_index][j_index][k_index] = u_hllc_single;
        } else {
          (void)u_hllc;
          auto const face_fluxes = riemann_solver_.SolveRiemannProblem<DIR>(
              material, reconstructed_conservatives_left,
              reconstructed_conservatives_right, reconstructed_primes_left,
              reconstructed_primes_right);
          for (unsigned int n = 0; n < MF::ANOE(); ++n) {
            fluxes[n][i_index][j_index][k_index] += face_fluxes[n];
          }
        }
      } // k
    }   // j
  }     // i
}
