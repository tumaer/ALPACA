//===---------------------- generic_stiffened_gas.h -----------------------===//
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
#ifndef GENERIC_STIFFENED_GAS_H
#define GENERIC_STIFFENED_GAS_H

#include "utilities/mathematical_functions.h"
#include <cmath>

namespace GenericStiffenedGas {
/**
 * @brief Computes pressure from inputs as -gamma*pi + (gamma - 1) * (E  - 0.5 *
 * rho * ||v^2||)
 * @tparam safe decision flag whether to take care of devision by 0
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Pressure according to stiffened-gas equation of state.
 */
template <bool safe>
constexpr double
CalculatePressure(const double density, const double momentum_x,
                  const double momentum_y, const double momentum_z,
                  const double energy, const double gamma, const double pi) {
  double const q = 0.5 * DimensionAwareConsistencyManagedSum(
                             momentum_x * momentum_x, momentum_y * momentum_y,
                             momentum_z * momentum_z);
  if constexpr (safe) {
    double const pressure =
        -gamma * pi +
        (gamma - 1.0) *
            (energy -
             q / std::max(density, std::numeric_limits<double>::epsilon()));
    return std::max(pressure, -pi + std::numeric_limits<double>::epsilon());
  }
  return -gamma * pi + (gamma - 1.0) * (energy - q / density);
}

/**
 * @brief Computes energy from inputs as (p + gamma * pi) / (gamma - 1) + 0.5 *
 * rho * ||v^2||
 * @param density .
 * @param velocity_x .
 * @param velocity_y .
 * @param velocity_z .
 * @param pressure .
 * @return Energy according to stiffened-gas equation of state.
 */
constexpr double CalculateEnergy(const double density, const double velocity_x,
                                 const double velocity_y,
                                 const double velocity_z, const double pressure,
                                 double const gamma, double const pi) {
  return (pressure + gamma * pi) / (gamma - 1.0) +
         (0.5 *
          DimensionAwareConsistencyManagedSum(velocity_x * velocity_x,
                                              velocity_y * velocity_y,
                                              velocity_z * velocity_z) *
          density);
}

/**
 * @brief Computes speed of sound from inputs as sqrt(gamma * (p + pi)) / rho
 * @tparam safe decision flag whether to take care of devision by 0 and negative
 * root
 * @param density .
 * @param pressure .
 * @return Speed of sound according to stiffened-gas equation of state.
 */
template <bool safe>
constexpr double CalculateSpeedOfSound(const double density,
                                       const double pressure,
                                       const double gamma, const double pi) {
  if constexpr (safe) {
    const double speed_of_sound_squared =
        gamma * (pressure + pi) /
        std::max(density, std::numeric_limits<double>::epsilon());
    return std::sqrt(std::max(speed_of_sound_squared,
                              std::numeric_limits<double>::epsilon()));
  } else {
    return std::sqrt(gamma * (pressure + pi) / density);
  }
}

} // namespace GenericStiffenedGas

#endif // GENERIC_STIFFENED_GAS_H
