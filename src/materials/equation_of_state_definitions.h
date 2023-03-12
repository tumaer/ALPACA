//===----------------- equation_of_state_definitions.h --------------------===//
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
#ifndef EQUATION_OF_STATE_DEFINITIONS_H
#define EQUATION_OF_STATE_DEFINITIONS_H

#include <stdexcept>
#include <string>

#include "utilities/string_operations.h"

/**
 * @brief The EquationOfStateName enum gives all EquationOfStates a unique
 * identifier to allow automatic selection of correct eos functions as well as
 * to distinguish treatment in case of contact between materials, e.g. merging
 * or interaction. No specific fixed index is required for this enum class,
 * since no indexing is done.
 */
enum class EquationOfStateName {
  StiffenedGas,
  StiffenedGasSafe,
  StiffenedGasCompleteSafe,
  NobleAbelStiffenedGas,
  WaterlikeFluid,
  Isentropic
};

/**
 * @brief Converts a string to its corresponding equation of state name.
 * @param eos_name The equation of state name.
 * @return Name of the equation of state.
 */
inline EquationOfStateName StringToEos(std::string const &eos_name) {
  // transform string to upper case without spaces
  std::string const eos_upper_case(
      StringOperations::ToUpperCaseWithoutSpaces(eos_name));
  // switch statements cannot be used with strings
  if (eos_upper_case == "STIFFENEDGAS") {
    return EquationOfStateName::StiffenedGas;
  } else if (eos_upper_case == "STIFFENEDGASSAFE") {
    return EquationOfStateName::StiffenedGasSafe;
  } else if (eos_upper_case == "STIFFENEDGASCOMPLETESAFE") {
    return EquationOfStateName::StiffenedGasCompleteSafe;
  } else if (eos_upper_case == "NOBLEABELSTIFFENEDGAS") {
    return EquationOfStateName::NobleAbelStiffenedGas;
  } else if (eos_upper_case == "WATERLIKEFLUID") {
    return EquationOfStateName::WaterlikeFluid;
  } else if (eos_upper_case == "ISENTROPIC") {
    return EquationOfStateName::Isentropic;
  } else {
    throw std::logic_error("Equation of state " + eos_upper_case +
                           " is not known!");
  }
}

#endif // EQUATION_OF_STATE_DEFINITIONS_H
