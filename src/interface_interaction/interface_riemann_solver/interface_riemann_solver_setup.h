//===------------------ interface_riemann_solver_setup.h ------------------===//
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
#ifndef INTERFACE_RIEMANN_SOLVER_SETUP_H
#define INTERFACE_RIEMANN_SOLVER_SETUP_H

#include "exact_iterative_interface_riemann_solver.h"
#include "hllc_interface_riemann_solver.h"
#include "linearized_interface_riemann_solver.h"
#include "two_rarefaction_iterative_interface_riemann_solver.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a InterfaceRiemannSolver type based on a specified
 * constexpr.
 */
namespace InterfaceRiemannSolverSetup {

/**
 * @brief Function returning the typedef of a InterfaceRiemannSolver based on a
 * constexpr template.
 *
 * @tparam InterfaceRiemannSolvers The constexpr template parameter to specify
 * the exact InterfaceRiemannSolver type.
 */
template <InterfaceRiemannSolvers> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<InterfaceRiemannSolvers::Linearized> {
  typedef LinearizedInterfaceRiemannSolver type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<InterfaceRiemannSolvers::Exact> {
  typedef ExactIterativeInterfaceRiemannSolver type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<InterfaceRiemannSolvers::TwoRarefaction> {
  typedef TwoRarefactionIterativeInterfaceRiemannSolver type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<InterfaceRiemannSolvers::Hllc> {
  typedef HllcInterfaceRiemannSolver type;
};

} // namespace InterfaceRiemannSolverSetup

#endif // INTERFACE_RIEMANN_SOLVER_SETUP_H
