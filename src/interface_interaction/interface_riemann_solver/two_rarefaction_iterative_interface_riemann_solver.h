//===-------- two_rarefaction_iterative_interface_riemann_solver.h --------===//
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
#ifndef TWO_RAREFACTION_ITERATIVE_INTERFACE_RIEMANN_SOLVER_H
#define TWO_RAREFACTION_ITERATIVE_INTERFACE_RIEMANN_SOLVER_H

#include "iterative_interface_riemann_solver.h"

/**
 * @brief Class for two rarefaction iterative interface riemann solver. Provides
 * functionality to compute the solution of a Riemann problem assuming a two
 * rarefaction solution.
 */
class TwoRarefactionIterativeInterfaceRiemannSolver
    : public IterativeInterfaceRiemannSolver<
          TwoRarefactionIterativeInterfaceRiemannSolver> {

  friend IterativeInterfaceRiemannSolver;

  std::array<double, 2> ObtainFunctionAndDerivativeImplementation(
      double const initial_root, double const p, double const pressure_function,
      double const one_pressure_function, double const pressure_constant,
      double const A, double const B, double const C, double const D) const;

public:
  TwoRarefactionIterativeInterfaceRiemannSolver() = delete;
  /**
   * @brief      Default constructor.
   *
   * @param[in]  material_manager  The material manager.
   */
  explicit TwoRarefactionIterativeInterfaceRiemannSolver(
      MaterialManager const &material_manager)
      : IterativeInterfaceRiemannSolver<
            TwoRarefactionIterativeInterfaceRiemannSolver>(
            material_manager) { /* Empty besides initialliser list */
  }
  ~TwoRarefactionIterativeInterfaceRiemannSolver() = default;
  TwoRarefactionIterativeInterfaceRiemannSolver(
      TwoRarefactionIterativeInterfaceRiemannSolver const &) = delete;
  TwoRarefactionIterativeInterfaceRiemannSolver &
  operator=(TwoRarefactionIterativeInterfaceRiemannSolver const &) = delete;
  TwoRarefactionIterativeInterfaceRiemannSolver(
      TwoRarefactionIterativeInterfaceRiemannSolver &&) = delete;
  TwoRarefactionIterativeInterfaceRiemannSolver &
  operator=(TwoRarefactionIterativeInterfaceRiemannSolver &&) = delete;
};

#endif // TWO_RAREFACTION_ITERATIVE_INTERFACE_RIEMANN_SOLVER_H
