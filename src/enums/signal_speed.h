//===-------------------------- signal_speed.h ----------------------------===//
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
#ifndef SIGNAL_SPEED_H
#define SIGNAL_SPEED_H

#include <string>

/**
 * @brief Identifier for nonlinear signal speed to be used in HLL-type Riemann
 * solvers.
 */
enum class SignalSpeed { Einfeldt, Davis, Toro, Arithmetic };

/**
 * @brief Provides a string representation of the applied SignalSpeed choice in
 * HLL-type Riemann solvers.
 * @param signal_speed SignalSpeed.
 * @return String representation for the given SignalSpeed choice.
 */
inline std::string SignalSpeedToString(SignalSpeed const signal_speed) {
  switch (signal_speed) {
  case SignalSpeed::Einfeldt:
    return "Einfeldt";
  case SignalSpeed::Davis:
    return "Davis";
  case SignalSpeed::Toro:
    return "Toro";
  case SignalSpeed::Arithmetic:
    return "Arithmetic";
  default:
    return "ERROR: This signal speed is not (yet) defined!";
  }
}

#endif // SIGNAL_SPEED_H
