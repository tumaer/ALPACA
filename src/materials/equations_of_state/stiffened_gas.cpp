//===------------------------ stiffened_gas.cpp ---------------------------===//
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
#include "materials/equations_of_state/stiffened_gas.h"
#include "materials/equations_of_state/generic_stiffened_gas.h"
#include "utilities/helper_functions.h"
#include "utilities/mathematical_functions.h"
#include "utilities/string_operations.h"
#include <cmath>

/**
 * @brief Constructs a stiffened gas equation of state with eos parameters given
 * as input.
 * @param dimensional_eos_data Map containing all data for the equation of
 * state.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 *
 * @note During the constructing a check is done if the required parameter
 * exists. If not an error is thrown. Furthermore, dimensionalization of each
 * value is done.
 */
StiffenedGas::StiffenedGas(
    std::unordered_map<std::string, double> const &dimensional_eos_data,
    UnitHandler const &unit_handler)
    : gamma_(GetCheckedParameter<double>(dimensional_eos_data, "gamma",
                                         "StiffenedGas")),
      background_pressure_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_eos_data,
                                      "backgroundPressure", "StiffenedGas"),
          UnitType::Pressure)) {
  /* Empty besides initializer list*/
}

/**
 * @brief Computes Pressure from inputs as -gamma*B + ( gamma - 1 ) * ( E  - 0.5
 * * rho * ||v^2|| ).
 * @param mass The mass used for the computation.
 * @param momentum_x The momentum in x-direction used for the computation.
 * @param momentum_y The momentum in y-direction used for the computation.
 * @param momentum_z The momentum in z-direction used for the computation.
 * @param energy The energy used for the computation.
 * @return Pressure according to stiffened-gas equation of state.
 */
double StiffenedGas::ComputePressure(double const mass, double const momentum_x,
                                     double const momentum_y,
                                     double const momentum_z,
                                     double const energy) const {
  return GenericStiffenedGas::CalculatePressure<false>(
      mass, momentum_x, momentum_y, momentum_z, energy, gamma_,
      background_pressure_);
}

/**
 * @brief Computes enthalpy as ( E + p ) / rho.
 * @param mass The mass used for the computation.
 * @param momentum_x The momentum in x-direction used for the computation.
 * @param momentum_y The momentum in y-direction used for the computation.
 * @param momentum_z The momentum in z-direction used for the computation.
 * @param energy The energy used for the computation.
 * @return Enthalpy value.
 */
double StiffenedGas::ComputeEnthalpy(double const mass, double const momentum_x,
                                     double const momentum_y,
                                     double const momentum_z,
                                     double const energy) const {
  return (energy +
          ComputePressure(mass, momentum_x, momentum_y, momentum_z, energy)) /
         mass;
}

/**
 * @brief Computes energy according to stiffened gas equation.
 * @param density The density used for the computation.
 * @param velocity_x The velocity in x-direction used for the computation.
 * @param velocity_y The velocity in y-direction used for the computation.
 * @param velocity_z The velocity in z-direction used for the computation.
 * @param pressure The pressure used for the computation.
 * @return Energy according to given inputs.
 */
double StiffenedGas::ComputeEnergy(double const density,
                                   double const velocity_x,
                                   double const velocity_y,
                                   double const velocity_z,
                                   double const pressure) const {
  return GenericStiffenedGas::CalculateEnergy(density, velocity_x, velocity_y,
                                              velocity_z, pressure, gamma_,
                                              background_pressure_);
}

/**
 * @brief Computes Gruneisen coefficient as ( gamma-1 ) for stiffened-gas
 * equation of state.
 * @return Gruneisen coefficient .
 */
double StiffenedGas::GetGruneisen() const { return (gamma_ - 1.0); }

/**
 * @brief Returns Gamma.
 * @return Gamma.
 */
double StiffenedGas::GetGamma() const { return gamma_; }

/**
 * @brief Returns B.
 * @return B.
 */
double StiffenedGas::GetB() const { return background_pressure_; }

/**
 * @brief Computes psi from inputs as ( p + gamma * B ) / rho.
 * @param pressure The pressure used for the computation.
 * @param one_density The density used for the computation.
 * @return Psi according to stiffened-gas equation of state.
 */
double StiffenedGas::ComputePsi(double const pressure,
                                double const one_density) const {
  return (pressure + gamma_ * background_pressure_) * one_density;
}

/**
 * @brief Computes Speed of Sound from inputs as sqrt( gamma * ( p + B ) ) /
 * rho.
 * @param density The density used for the computation.
 * @param pressure The pressure used for the computation.
 * @return Speed of sound according to stiffened-gas equation of state.
 */
double StiffenedGas::ComputeSpeedOfSound(double const density,
                                         double const pressure) const {
  return GenericStiffenedGas::CalculateSpeedOfSound<false>(
      density, pressure, gamma_, background_pressure_);
}

/**
 * @brief Provides logging information of the equation of state.
 * @param indent Number of white spaces used at the beginning of each line for
 * the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string StiffenedGas::GetLogData(unsigned int const indent,
                                     UnitHandler const &unit_handler) const {
  // string that is returned
  std::string log_string;
  // Name of the equation of state
  log_string += StringOperations::Indent(indent) +
                "Type                 : Stiffened gas\n";
  // Parameters with small indentation
  log_string +=
      StringOperations::Indent(indent) + "Gruneisen coefficient: " +
      StringOperations::ToScientificNotationString(GetGruneisen(), 9) + "\n";
  log_string += StringOperations::Indent(indent) + "Gamma                : " +
                StringOperations::ToScientificNotationString(gamma_, 9) + "\n";
  log_string += StringOperations::Indent(indent) + "Background pressure  : " +
                StringOperations::ToScientificNotationString(
                    unit_handler.DimensionalizeValue(background_pressure_,
                                                     UnitType::Pressure),
                    9) +
                "\n";
  return log_string;
}
