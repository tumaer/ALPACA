//===-------------- levelset_initializer_definitions.h --------------------===//
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
#ifndef LEVELSET_INITIALIZER_DEFINITIONS_H
#define LEVELSET_INITIALIZER_DEFINITIONS_H

#include <stdexcept>
#include <string>

#include "utilities/string_operations.h"

/**
 * @brief The LevelsetInitializerType enum gives all levelset initializer types
 * a unique identifier to allow automatic selection of correct levelset
 *        functions. No specific fixed index is required for this enum class,
 * since no indexing is done.
 */
enum class LevelsetInitializerType { Functional, STL, Parametric };

/**
 * @brief Converts a string to its corresponding levelset initializer type.
 * @param levelset_initializer_type The levelset initializer type.
 * @return Type of the levelset initializer.
 */
inline LevelsetInitializerType
StringToLevelsetInitializerType(std::string const &levelset_initializer_type) {
  // transform string to upper case without spaces
  std::string const levelset_initializer_upper_case(
      StringOperations::ToUpperCaseWithoutSpaces(levelset_initializer_type));
  // switch statements cannot be used with strings
  if (levelset_initializer_upper_case == "FUNCTIONAL") {
    return LevelsetInitializerType::Functional;
  } else if (levelset_initializer_upper_case == "STL") {
    return LevelsetInitializerType::STL;
  } else if (levelset_initializer_upper_case == "PARAMETRIC") {
    return LevelsetInitializerType::Parametric;
  } else {
    throw std::logic_error("Levelset initializer type " +
                           levelset_initializer_type + " is not known!");
  }
}

#endif // LEVELSET_INITIALIZER_DEFINITIONS_H
