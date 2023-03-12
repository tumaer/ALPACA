//===--------------------- dimension_definition.h -------------------------===//
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
#ifndef DIMENSION_DEFINITION_H
#define DIMENSION_DEFINITION_H

#include <type_traits>

/**
 * @brief Identifier for the dimensions to be simulated, i.e. one-dimensional,
 * two-dimensional or three-dimensional.
 *
 * @note It is absolutely neccesary to keep 1D = 1, 2D = 2 and 3D = 3.
 * Otherwise, the code will break. If additional Dimension types are required,
 * please append them. NEVER EVER change the underlying type.
 */
enum class Dimension : unsigned int { One = 1, Two = 2, Three = 3 };

/**
 * @brief Converts a dimension identifier to a (C++11 standard compliant, i. e.
 * positive) array index. "DTI = Dimension To Index".
 * @param d The dimension identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<Dimension>::type DTI(Dimension const d) {
  return static_cast<typename std::underlying_type<Dimension>::type>(d);
}

#endif // DIMENSION_DEFINITION_H
