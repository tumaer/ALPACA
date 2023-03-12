//===--------------- linearized_interface_riemann_solver.h ----------------===//
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
#ifndef LINEARIZED_INTERFACE_RIEMANN_SOLVER_H
#define LINEARIZED_INTERFACE_RIEMANN_SOLVER_H

#include "interface_riemann_solver.h"

/**
 * @brief Solves a Riemann problem linearized.
 */
class LinearizedInterfaceRiemannSolver
    : public InterfaceRiemannSolver<LinearizedInterfaceRiemannSolver> {

  friend InterfaceRiemannSolver;

  std::array<double, 3> SolveInterfaceRiemannProblemImplementation(
      double const rho_left, double const p_left,
      double const velocity_normal_left, MaterialName const material_left,
      double const rho_right, double const p_right,
      double const velocity_normal_right, MaterialName const material_right,
      double const delta_p) const;

public:
  LinearizedInterfaceRiemannSolver() = delete;
  explicit LinearizedInterfaceRiemannSolver(
      MaterialManager const &material_manager);
  LinearizedInterfaceRiemannSolver(LinearizedInterfaceRiemannSolver const &) =
      delete;
  LinearizedInterfaceRiemannSolver &
  operator=(LinearizedInterfaceRiemannSolver const &) = delete;
  LinearizedInterfaceRiemannSolver(LinearizedInterfaceRiemannSolver &&) =
      delete;
  LinearizedInterfaceRiemannSolver &
  operator=(LinearizedInterfaceRiemannSolver &&) = delete;
};

#endif // LINEARIZED_INTERFACE_RIEMANN_SOLVER_H
