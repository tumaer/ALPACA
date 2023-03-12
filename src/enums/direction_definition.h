//===--------------------- direction_definition.h -------------------------===//
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
#ifndef DIRECTION_DEFINITION_H
#define DIRECTION_DEFINITION_H

#include <type_traits>

/**
 * @brief Identifier for spatial direction, i.e. x, y and z.
 *
 * @note Do NOT change the indices of the enum entries (used for mapping in
 * arrays).
 */
enum class Direction : unsigned int { X = 0, Y = 1, Z = 2 };

/**
 * @brief Converts a direction identifier to a (C++11 standard compliant, i. e.
 * positive) array index. "DTI = Direction To Index".
 * @param d The direction identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<Direction>::type DTI(Direction const d) {
  return static_cast<typename std::underlying_type<Direction>::type>(d);
}

#endif // DIRECTION_DEFINITION_H
