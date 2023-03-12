//===----------------------- equation_settings.h --------------------------===//
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
#ifndef EQUATION_SETTINGS_H
#define EQUATION_SETTINGS_H

#include <string>

enum class EquationSet { Isentropic, Euler, NavierStokes, GammaModel, Custom };
constexpr EquationSet active_equations = EquationSet::NavierStokes;

enum class InterfaceSet { Default, Custom };
constexpr InterfaceSet active_interface_quantities = InterfaceSet::Default;

enum class ParameterSet { Default, Custom };
constexpr ParameterSet active_parameters = ParameterSet::Default;

/**
 * @brief provides a string representation of an equation set.
 * @param eq Equation set to be stringified.
 * @return String of the given equation set.
 */
inline std::string SetToString(EquationSet const eq) {
  switch (eq) {
  case EquationSet::Isentropic:
    return "Isentropic";
  case EquationSet::Euler:
    return "Euler";
  case EquationSet::NavierStokes:
    return "Navier-Stokes";
  case EquationSet::GammaModel:
    return "Gamma-Model";
  case EquationSet::Custom:
    return "Custom";
  default:
    return "ERROR: This set is not (yet) defined!";
  }
}

/**
 * @brief provides a string representation of an interface set.
 * @param eq Interface set to be stringified.
 * @return String of the given interface set.
 */
inline std::string SetToString(InterfaceSet const is) {
  switch (is) {
  case InterfaceSet::Default:
    return "Default";
  case InterfaceSet::Custom:
    return "Custom";
  default:
    return "ERROR: This set is not (yet) defined!";
  }
}

/**
 * @brief provides a string representation of a parameter set.
 * @param ps Parameter set to be stringified.
 * @return String of the given parameter set.
 */
inline std::string SetToString(ParameterSet const ps) {
  switch (ps) {
  case ParameterSet::Default:
    return "Default";
  case ParameterSet::Custom:
    return "Custom";
  default:
    return "ERROR: This set is not (yet) defined!";
  }
}

#endif // EQUATION_SETTINGS_H
