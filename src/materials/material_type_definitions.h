//===----------------- material_type_definitions.h ------------------------===//
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
#ifndef MATERIAL_TYPE_DEFINITIONS_H
#define MATERIAL_TYPE_DEFINITIONS_H

#include "utilities/string_operations.h"

/**
 * @brief The MaterialType enum gives all types a unique identifier.
 */
enum class MaterialType { Fluid, SolidBoundary };

/**
 * @brief Converts a string to its corresponding equation of state name.
 * @param eos_name The equation of state name.
 * @return Name of the equation of state.
 */
inline MaterialType
StringToMaterialType(std::string const &material_type_name) {
  // transform string to upper case without spaces
  std::string const type_upper_case(
      StringOperations::ToUpperCaseWithoutSpaces(material_type_name));
  // switch statements cannot be used with strings
  if (type_upper_case == "FLUID") {
    return MaterialType::Fluid;
  } else if (type_upper_case == "SOLIDBOUNDARY") {
    return MaterialType::SolidBoundary;
  } else {
    throw std::logic_error("Material type " + type_upper_case +
                           " is not known!");
  }
}

/**
 * @brief Gives the string for each material type.
 * @param type Material type for which the string should be returned.
 * @return string of material type.
 */
inline std::string MaterialTypeToString(MaterialType const type) {
  switch (type) {
  case MaterialType::Fluid: {
    return "Fluid";
  }
  case MaterialType::SolidBoundary: {
    return "Solid boundary";
  }
  default: {
    throw std::logic_error("Material type not known!");
  }
  }
}

#endif // MATERIAL_TYPE_DEFINITIONS_H
