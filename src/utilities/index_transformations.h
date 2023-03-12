//===-------------------- index_transformations.h -------------------------===//
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
#ifndef INDEX_TRANSFORMATIONS_H
#define INDEX_TRANSFORMATIONS_H

#include "user_specifications/compile_time_constants.h"

/**
 * @brief This namespace provides functionality to transform buffer indices.
 * "BIT" means "Buffer Index Transformation".
 */
namespace BIT {

/**
 * @brief Transforms an x-index in an buffer of size "total cells" to an index
 * in an buffer of size "internal cells". "T2I" means "Total Cell to Inner
 * Cell".
 * @param i The x-index in the buffer of size "total cells".
 * @return The x-index in the buffer of size "inner cells".
 */
constexpr unsigned int T2IX(unsigned int i) { return i - CC::FICX(); }

/**
 * @brief Transforms an y-index in an buffer of size "total cells" to an index
 * in an buffer of size "internal cells". "T2I" means "Total Cell to Inner
 * Cell".
 * @param j The y-index in the buffer of size "total cells".
 * @return The y-index in the buffer of size "inner cells".
 */
constexpr unsigned int T2IY(unsigned int j) {
  return CC::DIM() != Dimension::One ? j - CC::FICY() : 0;
}

/**
 * @brief Transforms an z-index in an buffer of size "total cells" to an index
 * in an buffer of size "internal cells". "T2I" means "Total Cell to Inner
 * Cell".
 * @param k The z-index in the buffer of size "total cells".
 * @return The z-index in the buffer of size "inner cells".
 */
constexpr unsigned int T2IZ(unsigned int k) {
  return CC::DIM() == Dimension::Three ? k - CC::FICZ() : 0;
}

/**
 * @brief Transforms an x-index in an buffer of size "inner cells" to an index
 * in an buffer of size "total cells". "I2T" means "Inner Cell to Total Cell".
 * @param i The x-index in the buffer of size "inner cells".
 * @return The x-index in the buffer of size "total cells".
 */
constexpr unsigned int I2TX(unsigned int i) { return i + CC::FICX(); }

/**
 * @brief Transforms an y-index in an buffer of size "inner cells" to an index
 * in an buffer of size "total cells". "I2T" means "Inner Cell to Total Cell".
 * @param j The y-index in the buffer of size "inner cells".
 * @return The y-index in the buffer of size "total cells".
 */
constexpr unsigned int I2TY(unsigned int j) {
  return CC::DIM() != Dimension::One ? j + CC::FICY() : 0;
}

/**
 * @brief Transforms an z-index in an buffer of size "inner cells" to an index
 * in an buffer of size "total cells". "I2T" means "Inner Cell to Total Cell".
 * @param k The z-index in the buffer of size "inner cells".
 * @return The z-index in the buffer of size "total cells".
 */
constexpr unsigned int I2TZ(unsigned int k) {
  return CC::DIM() == Dimension::Three ? k + CC::FICZ() : 0;
}

/**
 * @brief Transforms an x-index in an buffer of size "total cells" to an index
 * in an flux buffer. "T2F" means "Total Cell to Flux Array".
 * @param i The x-index in the buffer of size "total cells".
 * @return The x-index in the flux buffer.
 */
constexpr unsigned int T2FX(unsigned int i) { return i - (CC::FICX() - 1); }

/**
 * @brief Transforms an y-index in an buffer of size "total cells" to an index
 * in an flux buffer. "T2F" means "Total Cell to Flux Array".
 * @param j The y-index in the buffer of size "total cells".
 * @return The y-index in the flux buffer.
 */
constexpr unsigned int T2FY(unsigned int j) {
  return j - ((CC::DIM() != Dimension::One) ? (CC::FICY() - 1) : (-1));
}

/**
 * @brief Transforms an z-index in an buffer of size "total cells" to an index
 * in an flux buffer. "T2F" means "Total Cell to Flux Array".
 * @param k The z-index in the buffer of size "total cells".
 * @return The z-index in the flux buffer.
 */
constexpr unsigned int T2FZ(unsigned int k) {
  return k - ((CC::DIM() == Dimension::Three) ? (CC::FICZ() - 1) : (-1));
}

} // namespace BIT

#endif // INDEX_TRANSFORMATIONS_H
