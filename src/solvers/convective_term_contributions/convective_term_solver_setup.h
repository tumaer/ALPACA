//===------------------- convective_term_solver_setup.h -------------------===//
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
#ifndef CONVECTIVE_TERM_SOLVER_SETUP_H
#define CONVECTIVE_TERM_SOLVER_SETUP_H

#include "solvers/convective_term_contributions/finite_volume_scheme.h"
#include "solvers/convective_term_contributions/flux_splitting_scheme.h"
#include "user_specifications/riemann_solver_settings.h"

/**
 * @brief A namespace to get a ConvectiveTermSolver type based on a specified
 * constexpr.
 */
namespace ConvectiveTermSolverSetup {

/**
 * @brief Function returning the typedef of a ConvectiveTermSolvers based on a
 * constexpr template.
 *
 * @tparam ConvectiveTermSolvers The constexpr template parameter to specify the
 * exact ConvectiveTermSolvers type.
 */
template <ConvectiveTermSolvers> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ConvectiveTermSolvers::FluxSplitting> {
  typedef FluxSplittingScheme type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ConvectiveTermSolvers::FiniteVolume> {
  typedef FiniteVolumeScheme type;
};

} // namespace ConvectiveTermSolverSetup

#endif // CONVECTIVE_TERM_SOLVER_SETUP_H
