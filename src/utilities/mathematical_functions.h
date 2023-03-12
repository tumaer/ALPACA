//===--------------------- mathematical_functions.h -----------------------===//
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
#ifndef MATHEMATICAL_FUNCTIONS_H
#define MATHEMATICAL_FUNCTIONS_H

#include <algorithm>
#include <cmath>

#include "user_specifications/compile_time_constants.h"

/**
 * @brief Overload of signum function for positive only types ( e.g. unsigned
 * int ).
 */
template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, int>::type constexpr Signum(
    T const x) {
  return T(0) < x;
}

/**
 * @brief Computes the signum, i. e. -1,0,1 of the given input. Not to be
 * mistaken for sgn function which only gives +/- 1.
 */
template <typename T>
typename std::enable_if<std::is_signed<T>::value, int>::type constexpr Signum(
    T const x) {
  return (T(0) < x) - (x < T(0));
}

/**
 * @brief Computes the sign, i. e. -1,1 of the given input. Note, this functions
 * only gives +/- 1. Zero input results in 1 as result.
 */
template <typename T> T constexpr Sign(T const x) {
  return (T(0) <= x) - (x < T(0));
}

/**
 * @brief Computes the power function for the given input and an integer
 * exponent by repeated multiplication.
 */
template <unsigned int Exponent, typename T>
T constexpr IntegerPow([[maybe_unused]] T const x) {
  if constexpr (Exponent == 1) {
    return x;
  }
  if constexpr (Exponent > 1) {
    return IntegerPow<Exponent - 1>(x) * x;
  }
  // else: Exponent == 0
  return T(1);
}

/**
 * @brief Computes the floating point consistent sum of three double values. The
 * sum might be not the correctly rounded analytical one, however, the same
 * result is guaranteed when a, b and c are rearranged.
 */
constexpr double ConsistencyManagedSum(double const a, double const b,
                                       double const c) {
  if constexpr (CC::FUSY()) {
    double const sum1 = (a + b) + c;
    double const sum2 = (a + c) + b;
    double const sum3 = (b + c) + a;
    return 0.5 * (std::min({sum1, sum2, sum3}) + std::max({sum1, sum2, sum3}));
  } else {
    return a + b + c;
  }
}

/**
 * @brief Computes the floating point consistent sum of up to three double
 * values. Dependent on whether the simulation is 1D, 2D or 3D, only the first,
 * the first two or all three values are considered. For 3D: the sum might be
 * not the correctly rounded analytical one, however, the same result is
 * guaranteed when a, b and c are rearranged.
 */
constexpr double DimensionAwareConsistencyManagedSum(double const a,
                                                     double const b,
                                                     double const c) {
  if constexpr (CC::DIM() == Dimension::One) {
    return a;
  }
  if constexpr (CC::DIM() == Dimension::Two) {
    return a + b;
  }
  return ConsistencyManagedSum(a, b, c);
}

/**
 * @brief Computes the floating point consistent sum of four double values. The
 * sum might be not the correctly rounded analytical one, however, the same
 * result is guaranteed when a, b, c and d are rearranged.
 */
constexpr double ConsistencyManagedSum(double const a, double const b,
                                       double const c, double const d) {
  if constexpr (CC::FUSY()) {
    double const sum1 = (a + b) + (c + d);
    double const sum2 = (a + c) + (b + d);
    double const sum3 = (a + d) + (b + c);
    return 0.5 * (std::min({sum1, sum2, sum3}) + std::max({sum1, sum2, sum3}));
  } else {
    return a + b + c + d;
  }
}

/**
 * @brief Returns the value given in the array. Auxiliary function.
 */
constexpr double ConsistencyManagedSum(std::array<double, 1> const values) {
  return values[0];
}

/**
 * @brief Computes the floating point consistent sum of two double values in an
 * array. The sum might be not the correctly rounded analytical one, however,
 * the same result is guaranteed when a and b are rearranged.
 */
constexpr double ConsistencyManagedSum(std::array<double, 2> const values) {
  return values[0] + values[1];
}

/**
 * @brief Computes the floating point consistent sum of three double values in
 * an array. The sum might be not the correctly rounded analytical one, however,
 * the same result is guaranteed when a, b and c are rearranged.
 */
constexpr double ConsistencyManagedSum(std::array<double, 3> const values) {
  return ConsistencyManagedSum(values[0], values[1], values[2]);
}

/**
 * @brief Computes the floating point consistent sum of up to three double
 * values in an array. Dependent on whether the simulation is 1D, 2D or 3D, only
 * the first, the first two or all three values are considered. For 3D: the sum
 * might be not the correctly rounded analytical one, however, the same result
 * is guaranteed when a, b and c are rearranged.
 */
constexpr double
DimensionAwareConsistencyManagedSum(std::array<double, 3> const values) {
  return DimensionAwareConsistencyManagedSum(values[0], values[1], values[2]);
}

/**
 * @brief Computes the floating point consistent sum of four double values in an
 * array. The sum might be not the correctly rounded analytical one, however,
 * the same result is guaranteed when a, b, c and d are rearranged.
 */
constexpr double ConsistencyManagedSum(std::array<double, 4> const values) {
  return ConsistencyManagedSum(values[0], values[1], values[2], values[3]);
}

/**
 * @brief Implementation of the Min-Mod function according to \cite LeVeque1992.
 * @param value_1 The first value.
 * @param value_2 The second value.
 * @return The result of the Min-Mod operation.
 */
inline double MinMod(double const value_1, double const value_2) {
  return 0.5 * (Signum(value_1) + Signum(value_2)) *
         std::min(std::abs(value_1), std::abs(value_2));
}

/**
 * @brief Computes the absolute value.
 * @param value The value for which the absolute value should be obtained.
 * @return The absolute value.
 * @tparam T The type of the value.
 * @note Since the abs function of the STL library does not provide a constexpr
 * version, we implement our own constexpr noexcept version.
 */
template <typename T> constexpr T Abs(T const value) noexcept {
  return value >= static_cast<T>(0) ? value : (-value);
}

double GodunovHamiltonian(double const (&derivatives)[DTI(CC::DIM())][2],
                          double const old_levelset_sign);

#endif // MATHEMATICAL_FUNCTIONS_H
