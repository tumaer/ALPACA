//===----------------------- string_operations.h --------------------------===//
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
#ifndef STRING_OPERATIONS_H
#define STRING_OPERATIONS_H

#include <limits>
#include <string>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <type_traits>
#include <vector>

/**
 * @brief Provide the functionality to manipulate given strings
 */
namespace StringOperations {

// converts a number to scientific notation string, e.g. -1.0E+10
std::string ToScientificNotationString(
    double const number,
    int const precision = std::numeric_limits<double>::digits10 + 1,
    bool consider_sign = false);
// converts a string to a string without spaces and all letters are upper cases
std::string ToUpperCaseWithoutSpaces(std::string const &word);
// removes all spaces from a string
std::string RemoveSpaces(std::string const &word);
// gives a empty string with a certain width
std::string Indent(unsigned int const width);
// remove leading numbers from a string
std::string RemoveLeadingNumbers(std::string const &word);
// trim a string
std::string Trim(std::string const &word);

/**
 * @brief Converts a string into a variable of type T.
 * @param input The string that should be converted.
 * @tparam T Type of value to be returned.
 * @return The converted value.
 */
template <typename T> T ConvertStringToValue(std::string const input) {
  T return_value;
  std::stringstream temp(input);
  temp >> return_value;
  return return_value;
}

/**
 * @brief Converts string with white-space delimiter into vector of type T.
 * @param input_data String to be converted to vector.
 * @tparam T Type of vector to be returned.
 * @return Vector of values.
 */
template <typename T>
std::vector<T> ConvertStringToVector(std::string const input_data) {

  // cut string at white spaces
  std::istringstream string_iterator(input_data);
  std::vector<std::string> const string_vector(
      std::istream_iterator<std::string>{string_iterator},
      std::istream_iterator<std::string>());

  // return vector if type is string, otherwise convert to requested type and
  // then return
  if constexpr (std::is_same<T, std::string>::value) {
    return string_vector;
  } else {
    std::vector<T> return_vector(string_vector.size());
    std::transform(std::cbegin(string_vector), std::cend(string_vector),
                   std::begin(return_vector), [](std::string const &str) {
                     return ConvertStringToValue<T>(str);
                   });
    return return_vector;
  }
}
} // namespace StringOperations

#endif // STRING_OPERATIONS_H
