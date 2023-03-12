//===------------------------- riemann_solver.h ---------------------------===//
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
#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H
#include <array>

#include "block_definitions/block.h"
#include "materials/material_manager.h"
#include "solvers/eigendecomposition.h"

/**
 * @brief Interface to solve the underlying system of equations. Uses spatial
 * reconstruction stencils to approximate the solution.
 */
template <typename DerivedRiemannSolver> class RiemannSolver {

  friend DerivedRiemannSolver;

  EigenDecomposition const &eigendecomposition_calculator_;
  MaterialManager const &material_manager_;

  explicit RiemannSolver(
      MaterialManager const &material_manager,
      EigenDecomposition const &eigendecomposition_calculator)
      : eigendecomposition_calculator_(eigendecomposition_calculator),
        material_manager_(material_manager) {
    // Empty constructor besides initializer list.
  }

public:
  RiemannSolver() = delete;
  ~RiemannSolver() = default;
  RiemannSolver(RiemannSolver const &) = delete;
  RiemannSolver &operator=(RiemannSolver const &) = delete;
  RiemannSolver(RiemannSolver &&) = delete;
  RiemannSolver &operator=(RiemannSolver &&) = delete;

  /**
   * @brief Solves the first-order Riemann problem using left and right state
   * vectors.
   * @tparam direction.
   * @param material Container holding the relevant fluid data to compute the
   * update.
   * @param state_face_left, state_face_right The reconstructed left/right
   * conservative states.
   * @param prime_state_left, prime_state_right The primitive left/right
   * conservative states.
   * @return The flux over the cell face as computed by this Riemann solver.
   */
  template <Direction DIR>
  std::array<double, MF::ANOE()> SolveRiemannProblem(
      MaterialName const material,
      std::array<double, MF::ANOE()> const &state_face_left,
      std::array<double, MF::ANOE()> const &state_face_right,
      std::array<double, MF::ANOP()> const &prime_state_left,
      std::array<double, MF::ANOP()> const &prime_state_right) const {
    return static_cast<DerivedRiemannSolver const &>(*this)
        .template SolveRiemannProblemImplementation<DIR>(
            material, state_face_left, state_face_right, prime_state_left,
            prime_state_right);
  }

  /**
   * @brief Solves the first-order Riemann problem using left and right state
   * vectors.
   * @tparam direction.
   * @param material Container holding the relevant fluid data to compute the
   * update.
   * @param state_face_left, state_face_right The reconstructed left/right
   * conservative states.
   * @param prime_state_left, prime_state_right The primitive left/right
   * conservative states.
   * @return The flux over the cell face as computed by this Riemann solver and
   * the reconstructed velocity u_hllc (required for source terms).
   */
  template <Direction DIR>
  std::tuple<std::array<double, MF::ANOE()>, double> SolveGammaRiemannProblem(
      MaterialName const material,
      std::array<double, MF::ANOE()> const &state_face_left,
      std::array<double, MF::ANOE()> const &state_face_right,
      std::array<double, MF::ANOP()> const &prime_state_left,
      std::array<double, MF::ANOP()> const &prime_state_right) const {
    return static_cast<DerivedRiemannSolver const &>(*this)
        .template SolveGammaRiemannProblemImplementation<DIR>(
            material, state_face_left, state_face_right, prime_state_left,
            prime_state_right);
  }
};

#endif // RIEMANN_SOLVER_H
