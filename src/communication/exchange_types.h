//===------------------------ exchange_types.h ----------------------------===//
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
#ifndef EXCHANGE_TYPES_H
#define EXCHANGE_TYPES_H

#include "user_specifications/compile_time_constants.h"

/**
 * @brief Enum class to give the index of the appropriate Exchanges Types.
 * @note Do not change the underlying type. Used for index mapping.
 */
enum class ExchangeType : unsigned short { Plane = 0, Stick = 1, Cube = 2 };

/**
 * @brief Converts an ExchangeType identifier to a (C++11 standard compliant, i.
 * e. positive) array index. "ETTI = Exchange Type To Index".
 * @param et The exchange type identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<ExchangeType>::type ETTI(ExchangeType const et) {
  return static_cast<typename std::underlying_type<ExchangeType>::type>(et);
}

/**
 * @brief Exchange type for a single plane (east, west, north, south, top,
 * bottom).
 */
struct ExchangePlane {
  double plane_[CC::HSSX()][CC::ICY()][CC::ICZ()];
};
static_assert(sizeof(ExchangePlane) ==
                  CC::HSSX() * CC::ICY() * CC::ICZ() * sizeof(double),
              "ExchangePlane Struct is not contiguous in Memory");

/**
 * @brief Exchange type for a single stick representing the intersection between
 * two planes (e.g., north-east).
 */
struct ExchangeStick {
  double stick_[CC::HSSX()][CC::HSSY()][CC::ICZ()];
};
static_assert(sizeof(ExchangeStick) ==
                  CC::HSSX() * CC::HSSY() * CC::ICZ() * sizeof(double),
              "ExchangeStick Struct is not contiguous in Memory");

/**
 * @brief Exchange type for a single cube representing the intersection between
 * three planes (e.g., north-east-bottom).
 */
struct ExchangeCube {
  double cube_[CC::HSSX()][CC::HSSY()][CC::HSSZ()];
};
static_assert(sizeof(ExchangeCube) ==
                  CC::HSSX() * CC::HSSY() * CC::HSSZ() * sizeof(double),
              "ExchangeCube Struct is not contiguous in Memory");

/**
 * @brief Container of the interface tag array. In order to store it in
 * std::vectors. Used for more efficient MPI communication.
 */
struct InterfaceTagBundle {
  std::int8_t interface_tags_[CC::TCX()][CC::TCY()][CC::TCZ()];
};
// Check Memory Layout at compile time for safe MPI sending (Ensures Compiler
// did not pad the struct)
static_assert(sizeof(InterfaceTagBundle) ==
                  CC::TCX() * CC::TCY() * CC::TCZ() * sizeof(std::int8_t),
              "InterfaceTagBundle is not contiguous in Memory");

#endif // EXCHANGE_TYPES_H
