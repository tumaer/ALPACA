//===---------------------- input_definitions.h ---------------------------===//
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
#ifndef INPUT_DEFINITIONS_H
#define INPUT_DEFINITIONS_H

#include "utilities/string_operations.h"
#include <stdexcept>
#include <string>

/**
 * @brief The InputType enum defines the type of the input file that is used for
 * the output. (XML: An xml file is used as input ).
 */
enum class InputType { Xml };

/**
 * @brief Converts the InputType to its corresponding string (for logging).
 * @param type The Input type identifier.
 * @return String to be used.
 */
inline std::string InputTypeToString(InputType const type) {

  switch (type) {
  case InputType::Xml: {
    return "Xml";
  }
  default: {
    throw std::logic_error("Input type is not known!");
  }
  }
}

/**
 * @brief Gives the input type for a given string.
 * @param type Type as a string.
 * @return input type identifier.
 */
inline InputType StringToInputType(std::string const &type) {
  // Get the file extension
  std::string const type_upper_case(
      StringOperations::ToUpperCaseWithoutSpaces(type));
  // cannot use switch for strings
  if (type_upper_case == "XML") {
    return InputType::Xml;
  } else {
    // throw error if it is does not exist
    throw std::logic_error("Input file type '" + type_upper_case +
                           "' is not known!");
  }
}

#endif // INPUT_DEFINITIONS_H
