//===--------------------- string_operations.cpp --------------------------===//
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
#include "utilities/string_operations.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

namespace StringOperations {

/**
 * @brief Gives a string representation of the given input in scientific
 * notation.
 * @param number Floating-point value.
 * @param precision Allows to control the number of displayed decimals.
 * @param consider_sign Flag if an extra space should be considered for a
 * possible negative sign of positive numbers.
 * @return String-converted value.
 */
std::string ToScientificNotationString(double const number, int const precision,
                                       bool consider_sign) {
  // If number is positive add extra space for negative sign
  std::string sign_string = "";
  if (consider_sign && number >= 0.0) {
    sign_string = " ";
  }
  std::ostringstream out;
  out << sign_string << std::scientific << std::setprecision(precision)
      << number;
  return out.str();
}

/**
 * @brief Converts a string into an upper-case string without white spaces.
 * @param word String that has to be converted.
 * @return converted string.
 */
std::string ToUpperCaseWithoutSpaces(std::string const &word) {
  // Make it local due to erasing
  std::string upper_case_word(word);
  // Remove white spaces
  upper_case_word.erase(
      std::remove_if(upper_case_word.begin(), upper_case_word.end(), ::isspace),
      upper_case_word.end());
  // Convert to upper case
  std::transform(upper_case_word.begin(), upper_case_word.end(),
                 upper_case_word.begin(), ::toupper);

  return upper_case_word;
}

/**
 * @brief Converts a string into a string without white spaces.
 * @param word String that has to be converted.
 * @return converted string.
 */
std::string RemoveSpaces(std::string const &word) {
  // Make it local due to erasing
  std::string word_without_spaces(word);
  // Remove white spaces
  word_without_spaces.erase(std::remove_if(word_without_spaces.begin(),
                                           word_without_spaces.end(),
                                           ::isspace),
                            word_without_spaces.end());

  return word_without_spaces;
}

/**
 * @brief Removes all leading and trailing whitespaces, newline characters and
 * tabs from a string.
 * @param word String that has to be converted.
 * @return trimmed string.
 */
std::string Trim(std::string const &word) {
  // Make it local due to erasing
  std::string word_without_spaces(word);
  // Get the first element that does not contain any of the deilimiters
  auto const substring_start = word.find_first_not_of(" \n\t");
  if (substring_start == std::string::npos) {
    return "";
  }
  auto const substring_end = word.find_last_not_of(" \n\t");
  return word.substr(substring_start, substring_end - substring_start + 1);
}

/**
 * @brief Gives an empty string of certain width.
 * @param width Width of the string.
 * @return String.
 */
std::string Indent(unsigned int const width) { return std::string(width, ' '); }

/**
 * @brief Removes the leading numbers of as string.
 * @param word The string that should be modified.
 * @return The modified string.
 */
std::string RemoveLeadingNumbers(std::string const &word) {
  std::string name_without_leading_number(word);
  return name_without_leading_number.erase(
      0, std::min(word.find_first_not_of("0123456789"), word.size() - 1));
}
} // namespace StringOperations
