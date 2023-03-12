//===----------------------- vector_utilities.h ---------------------------===//
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
#ifndef VECTOR_UTILITIES_H
#define VECTOR_UTILITIES_H

#include <algorithm>
#include <functional>
#include <numeric>

#include "utilities/mathematical_functions.h"

namespace VectorUtilities {

/**
 * @brief Multiplies a single factor on all elements of the given vector.
 * @param vector The vector that should be multiplied.
 * @param factor The factor to be multiplied.
 * @return The mutliplied vector.
 */
inline std::array<double, 3> Multiply(std::array<double, 3> const vector,
                                      double const factor) {
  std::array<double, 3> return_array = {0.0, 0.0, 0.0};
  std::transform(
      std::cbegin(vector), std::cend(vector), std::begin(return_array),
      std::bind(std::multiplies<double>(), std::placeholders::_1, factor));
  return return_array;
}

/**
 * @brief Computes the L2 norm of a vector.
 * @param vector Vector for which the L2 norm should be computed.
 * @return L2 norm of the vector.
 * @note The computation consideres only the active dimensions.
 */
inline double L2Norm(std::array<double, 3> const vector) {
  return std::sqrt(DimensionAwareConsistencyManagedSum(
      vector[0] * vector[0], vector[1] * vector[1], vector[2] * vector[2]));
}

/**
 * @brief Computes the cross product between two vectors.
 * @param vector1, vector2 Two vectors for which the cross product should be
 * computed.
 * @return cross product between the two vectors.
 * @note Immutable operation. A flip of both vectors results in the same value,
 * but opposite sign.
 */
constexpr std::array<double, 3>
CrossProduct(std::array<double, 3> const vector1,
             std::array<double, 3> const vector2) {
  return {vector1[1] * vector2[2] - vector1[2] * vector2[1],
          vector1[2] * vector2[0] - vector1[0] * vector2[2],
          vector1[0] * vector2[1] - vector1[1] * vector2[0]};
}

/**
 * @brief Computes the dot product of two vectors.
 * @param vector1, vector2 Two vectors for which the dot product should be
 * computed.
 * @return dot product of the two vectors.
 * @note Mutable operation.
 */
inline double DotProduct(std::array<double, 3> const vector1,
                         std::array<double, 3> const vector2) {
  return std::inner_product(std::cbegin(vector1), std::cend(vector1),
                            std::cbegin(vector2), 0.0);
}

/**
 * @brief Computes the scalar triple product of three vectors in the order:
 * (vector1 x vector2) vector3
 * @param vector1, vector2, vector3 The three vectors of the triple product.
 * @return The scalar triple product value.
 * @note A circular shift of the vectors gives the same result. Swapping two
 * vectors leads to the opposite sign.
 */
inline double ScalarTripleProduct(std::array<double, 3> const vector1,
                                  std::array<double, 3> const vector2,
                                  std::array<double, 3> const vector3) {
  return DotProduct(CrossProduct(vector1, vector2), vector3);
}

/**
 * @brief Computes the vectorial triple product of three vectors in the order:
 * (vector1 x vector2) x vector3
 * @param vector1, vector2, vector3 The three vectors of the triple product.
 * @return The vectorial triple product value.
 * @note A swap of the first two vectors leads to a sign change.
 */
constexpr std::array<double, 3>
VectorTripleProduct(std::array<double, 3> const vector1,
                    std::array<double, 3> const vector2,
                    std::array<double, 3> const vector3) {
  return CrossProduct(CrossProduct(vector1, vector2), vector3);
}

/**
 * @brief Normalizes a vector.
 * @param input vector which is normalized.
 * @return The normalized vector.
 */
inline std::array<double, 3> Normalize(std::array<double, 3> const vector) {
  return Multiply(vector, 1.0 / L2Norm(vector));
}

/**
 * @brief Computes difference of two vectors (vector2 - vector1).
 * @param vector1, vector2 Vectors which are subtracted.
 * @return Difference vector.
 */
inline std::array<double, 3> Difference(std::array<double, 3> const vector1,
                                        std::array<double, 3> const vector2) {
  std::array<double, 3> return_array = {0.0, 0.0, 0.0};
  std::transform(std::cbegin(vector2), std::cend(vector2), std::cbegin(vector1),
                 std::begin(return_array), std::minus<double>());
  return return_array;
}

/**
 * @brief Computes distance between two vectors.
 * @param vector1, vector2 Vectors for which the distance is computed.
 * @return The distance.
 * @note Mutable operation.
 */
inline double Distance(std::array<double, 3> const vector1,
                       std::array<double, 3> const vector2) {
  return L2Norm(Difference(vector1, vector2));
}

} // namespace VectorUtilities

namespace VU = VectorUtilities;

#endif // VECTOR_UTILITIES_H
