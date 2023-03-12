//===------------ exact_iterative_interface_riemann_solver.cpp ------------===//
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
#include "exact_iterative_interface_riemann_solver.h"

/**
 * @brief      Calculates the shock/rarefaction relation and its derivative. See
 * \cite Toro2009.
 *
 * @param[in]  initial_root           The interface pressure of the current
 * iteration.
 * @param[in]  p                      The left/right pressure defining the
 * interface Riemann problem.
 * @param[in]  pressure_function      A constant dependent on the EoS.
 * @param[in]  one_pressure_function  The inverse of the pressure function.
 * @param[in]  pressure_constant      The pressure constant B.
 * @param[in]  A                      A constant. See IterationConstants::A()
 * for details.
 * @param[in]  B                      A constant. See IterationConstants::B()
 * for details.
 * @param[in]  C                      A constant. See IterationConstants::C()
 * for details.
 * @param[in]  D                      A constant. See IterationConstants::D()
 * for details.
 *
 * @return     The shock/rarefaction relation and its derivative.
 */
std::array<double, 2>
ExactIterativeInterfaceRiemannSolver::ObtainFunctionAndDerivativeImplementation(
    double const initial_root, double const p, double const pressure_function,
    double const one_pressure_function, double const pressure_constant,
    double const A, double const B, double const C, double const D) const {

  if (initial_root > p) {
    return {
        ShockRelations::Function(IterationUtilities::MaterialPressureFunction(
                                     initial_root, pressure_constant),
                                 pressure_function, A, B),
        ShockRelations::Derivative(IterationUtilities::MaterialPressureFunction(
                                       initial_root, pressure_constant),
                                   pressure_function, A, B)};
  } else {
    return {RarefactionRelations::Function(
                IterationUtilities::MaterialPressureFunction(initial_root,
                                                             pressure_constant),
                one_pressure_function, C, D),
            RarefactionRelations::Derivative(
                IterationUtilities::MaterialPressureFunction(initial_root,
                                                             pressure_constant),
                one_pressure_function, C, D)};
  }
}
