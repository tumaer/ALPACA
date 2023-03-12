//===----------------- hllc_interface_riemann_solver.h --------------------===//
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
#ifndef HLLC_INTERFACE_RIEMANN_SOLVER_H
#define HLLC_INTERFACE_RIEMANN_SOLVER_H

#include "interface_riemann_solver.h"

/**
 * @brief Solves a Riemann problem using a generalized HLLC approach according
 * to \cite Hu2008.
 * @note NOT USABLE TOGETHER WITH WATERLIKE EOS.
 * @note NOT USABLE FOR SURFACE TENSION COMPUTATIONS.
 */
class HllcInterfaceRiemannSolver
    : public InterfaceRiemannSolver<HllcInterfaceRiemannSolver> {

  friend InterfaceRiemannSolver;

  std::array<double, 3> SolveInterfaceRiemannProblemImplementation(
      double const rho_left, double const p_left,
      double const velocity_normal_left, MaterialName const material_left,
      double const rho_right, double const p_right,
      double const velocity_normal_right, MaterialName const material_right,
      double const delta_p) const;

public:
  explicit HllcInterfaceRiemannSolver() = delete;
  HllcInterfaceRiemannSolver(MaterialManager const &material_manager);
  ~HllcInterfaceRiemannSolver() = default;
  HllcInterfaceRiemannSolver(HllcInterfaceRiemannSolver const &) = delete;
  HllcInterfaceRiemannSolver &
  operator=(HllcInterfaceRiemannSolver const &) = delete;
  HllcInterfaceRiemannSolver(HllcInterfaceRiemannSolver &&) = delete;
  HllcInterfaceRiemannSolver &operator=(HllcInterfaceRiemannSolver &&) = delete;
};

#endif // HLLC_INTERFACE_RIEMANN_SOLVER_H
