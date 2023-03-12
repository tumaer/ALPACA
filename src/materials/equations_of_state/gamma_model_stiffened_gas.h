//===------------------- gamma_model_stiffened_gas.h ----------------------===//
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
#ifndef GAMMA_MODEL_STIFFENED_H
#define GAMMA_MODEL_STIFFENED_H

#include "materials/equations_of_state/generic_stiffened_gas.h"
#include "user_specifications/gamma_model_settings.h"
#include "utilities/mathematical_functions.h"
#include <cmath>

namespace GammaModelStiffenedGas {

/**
 * @brief Computes pressure from inputs as -gamma*pi + (gamma - 1) * (E  - 0.5 *
 * rho * ||v^2||)
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Pressure according to stiffened-gas equation of state.
 */
constexpr double CalculatePressure(const double density,
                                   const double momentum_x,
                                   const double momentum_y,
                                   const double momentum_z, const double energy,
                                   const double gamma, const double pi) {
  return GenericStiffenedGas::CalculatePressure<
      GammaModelSettings::EosSafeGuarding>(density, momentum_x, momentum_y,
                                           momentum_z, energy, gamma, pi);
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
  return GenericStiffenedGas::CalculateEnergy(density, velocity_x, velocity_y,
                                              velocity_z, pressure, gamma, pi);
}

/**
 * @brief Computes speed of sound from inputs as sqrt(gamma * (p + pi)) / rho
 * @param density .
 * @param pressure .
 * @return Speed of sound according to stiffened-gas equation of state.
 */
constexpr double CalculateSpeedOfSound(const double density,
                                       const double pressure,
                                       const double gamma, const double pi) {
  return GenericStiffenedGas::CalculateSpeedOfSound<
      GammaModelSettings::EosSafeGuarding>(density, pressure, gamma, pi);
}

/**
 * @brief Computes Gamma from inputs as 1 / (gamma - 1)
 * @param gamma .
 * @return Gamma according to Gamma Model requirements.
 */
constexpr double CalculateGamma(double const gamma) {
  return 1.0 / (gamma - 1.0);
}

/**
 * @brief Computes gamma from inputs as 1 / Gamma + 1
 * @param Gamma .
 * @return gamma according to Gamma Model requirements.
 */
constexpr double CalculatePrimeGamma(double const Gamma) {
  return 1.0 / Gamma + 1.0;
}

/**
 * @brief Computes Pi from inputs as ( gamma / ( gamma - 1 ) ) * pi
 * @param gamma .
 * @param pi .
 * @return Pi according to Gamma Model requirements.
 */
constexpr double CalculatePi(double const gamma, double const pi) {
  return (gamma / (gamma - 1.0)) * pi;
}

/**
 * @brief Computes pi from inputs as ( ( gamma - 1 ) / gamma ) * Pi
 * @param gamma .
 * @param Pi .
 * @return pi according to Gamma Model requirements.
 */
constexpr double CalculatePrimePi(double const gamma, double const Pi) {
  return ((gamma - 1.0) / gamma) * Pi;
}
} // namespace GammaModelStiffenedGas

#endif // GAMMA_MODEL_STIFFENED_H
