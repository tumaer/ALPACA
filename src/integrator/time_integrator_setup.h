//===--------------------- time_integrator_setup.h ------------------------===//
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
#ifndef TIME_INTEGRATOR_SETUP_H
#define TIME_INTEGRATOR_SETUP_H

#include "runge_kutta_2_TVD.h"
#include "runge_kutta_3_TVD.h"
#include "time_integrator.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a TimeIntegrator type based on a specified
 * constexpr.
 */
namespace TimeIntegratorSetup {

/**
 * @brief Function returning the typedef of a TimeIntegrator based on a
 * constexpr template.
 *
 * @tparam TimeIntegrators The constexpr template parameter to specify the exact
 * TimeIntegrator type.
 */
template <TimeIntegrators> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<TimeIntegrators::RK2> {
  typedef RungeKutta2TVD type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<TimeIntegrators::RK3> {
  typedef RungeKutta3TVD type;
};

} // namespace TimeIntegratorSetup

#endif // TIME_INTEGRATOR_SETUP_H
