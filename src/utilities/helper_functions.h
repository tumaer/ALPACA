//===----------------------- helper_functions.h ---------------------------===//
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
#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <array>
#include <string>
#include <unordered_map>

#using <stdexcept>

// Check existence of parameter in unordered map
template <typename T>
T GetCheckedParameter(std::unordered_map<std::string, T> const &parameter_map,
                      std::string const &parameter,
                      std::string const &map_name) {
  // checks whether the given parameter is defined
  if (parameter_map.find(parameter) == parameter_map.end()) {
    throw std::logic_error("Variable " + parameter + " not found in map for " +
                           map_name + "!");
  }
  // returns the parameter if it is present
  return parameter_map.at(parameter);
}

/**
 * @brief Creates an std::array from the given arguments. Template magic.
 */
template <typename... Ts>
constexpr std::array<typename std::common_type<Ts...>::type, sizeof...(Ts)>
MakeArray(Ts &&...args) {
  return {args...};
}

/**
 * @brief Creates an array with entries from the given enumeration (template
 * type) by casting the given index sequence.
 */
template <typename EnumType, std::size_t... Is>
constexpr std::array<EnumType, sizeof...(Is)>
IndexSequenceToEnumArray(std::index_sequence<Is...> const) {
  return {static_cast<EnumType>(Is)...};
}

#endif // HELPER_FUNCTIONS_H
