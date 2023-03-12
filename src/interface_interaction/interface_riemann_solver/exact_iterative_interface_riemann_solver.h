//===---------- exact_iterative_interface_riemann_solver.h ----------------===//
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
#ifndef EXACT_ITERATIVE_INTERFACE_RIEMANN_SOLVER_H
#define EXACT_ITERATIVE_INTERFACE_RIEMANN_SOLVER_H

#include "iterative_interface_riemann_solver.h"

/**
 * @brief Class for exact iterative interface riemann solver. Computes the exact
 * solution of a Riemann problem.
 */
class ExactIterativeInterfaceRiemannSolver
    : public IterativeInterfaceRiemannSolver<
          ExactIterativeInterfaceRiemannSolver> {

  friend IterativeInterfaceRiemannSolver;

  std::array<double, 2> ObtainFunctionAndDerivativeImplementation(
      double const initial_root, double const p, double const pressure_function,
      double const one_pressure_function, double const pressure_constant,
      double const A, double const B, double const C, double const D) const;

public:
  ExactIterativeInterfaceRiemannSolver() = delete;
  /**
   * @brief Default constructor.
   *
   * @param[in]  material_manager  The material manager.
   */
  ExactIterativeInterfaceRiemannSolver(MaterialManager const &material_manager)
      : IterativeInterfaceRiemannSolver<ExactIterativeInterfaceRiemannSolver>(
            material_manager) {}
  ~ExactIterativeInterfaceRiemannSolver() = default;
  ExactIterativeInterfaceRiemannSolver(
      ExactIterativeInterfaceRiemannSolver const &) = delete;
  ExactIterativeInterfaceRiemannSolver
  operator=(ExactIterativeInterfaceRiemannSolver const &) = delete;
  ExactIterativeInterfaceRiemannSolver(
      ExactIterativeInterfaceRiemannSolver &&) = delete;
  ExactIterativeInterfaceRiemannSolver &
  operator=(ExactIterativeInterfaceRiemannSolver &&) = delete;
};

#endif // EXACT_ITERATIVE_INTERFACE_RIEMANN_SOLVER_H
