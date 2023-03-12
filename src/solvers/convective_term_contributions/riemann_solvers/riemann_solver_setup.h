//===----------------------- riemann_solver_setup.h -----------------------===//
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
#ifndef RIEMANN_SOLVER_SETUP_H
#define RIEMANN_SOLVER_SETUP_H

#include "solvers/convective_term_contributions/riemann_solvers/gamma_hllc_lm_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/gamma_hllc_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/hll_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/hllc_lm_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/hllc_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/isentropic_hll_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/isentropic_hllc_riemann_solver.h"
#include "user_specifications/equation_settings.h"
#include "user_specifications/numerical_setup.h"
#include "user_specifications/riemann_solver_settings.h"

/**
 * @brief A namespace to get a RiemannSolver type based on a specified
 * constexpr.
 */
namespace RiemannSolverSetup {

namespace Isentropic {
/**
 * @brief Function returning the Isentropic Riemann solver matching the type in
 * the template argument.
 * @tparam RiemannSolvers Specification of the RiemannSolver type.
 */
template <FiniteVolumeSettings::RiemannSolvers> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc> {
  using type = IsentropicHllcRiemannSolver;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hll> {
  using type = IsentropicHllRiemannSolver;
};
} // namespace Isentropic

namespace EulerNavierStokes {
/**
 * @brief Function returning the Euler or Navier-Stokes equations Riemann solver
 * matching the type in the template argument.
 * @tparam RiemannSolvers Specification of the RiemannSolver type.
 */
template <FiniteVolumeSettings::RiemannSolvers> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc> {
  using type = HllcRiemannSolver;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc_LM> {
  using type = HllcLMRiemannSolver;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hll> {
  using type = HllRiemannSolver;
};
} // namespace EulerNavierStokes

namespace GammaModel {
/**
 * @brief Function returning the Gamma Model Riemann solver matching the type in
 * the template argument.
 * @tparam RiemannSolvers Specification of the RiemannSolver type.
 */
template <FiniteVolumeSettings::RiemannSolvers> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc> {
  using type = GammaHllcRiemannSolver;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc_LM> {
  using type = GammaHllcLMRiemannSolver;
};
} // namespace GammaModel

/**
 * @brief Function returning the Riemann solver matching the type in the
 * template argument accroding to the equation set in the second template
 * parameter.
 * @tparam RiemannSolvers Specification of the RiemannSolver type.
 * @tparam EquationSet Specification of the Equation(s) beeing solved.
 */
template <FiniteVolumeSettings::RiemannSolvers R, EquationSet> struct Dispatch {
  using type = typename EulerNavierStokes::Concretize<R>::type;
};

/**
 * @brief See generic implementation.
 */
template <FiniteVolumeSettings::RiemannSolvers R>
struct Dispatch<R, EquationSet::Isentropic> {
  using type = typename Isentropic::Concretize<R>::type;
};

/**
 * @brief See generic implementation.
 */
template <FiniteVolumeSettings::RiemannSolvers R>
struct Dispatch<R, EquationSet::GammaModel> {
  using type = typename GammaModel::Concretize<R>::type;
};

/**
 * @brief Function returning the Riemann solver for the (globally) selected
 * Equation and the given Solver template argument.
 * @tparam RiemannSolvers Specification of the RiemannSolver type.
 */
template <FiniteVolumeSettings::RiemannSolvers R> struct Concretize {
  using type = typename Dispatch<R, active_equations>::type;
};

} // namespace RiemannSolverSetup

#endif // RIEMANN_SOLVER_SETUP_H
