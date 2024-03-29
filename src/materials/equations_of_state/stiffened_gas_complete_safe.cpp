//===------------------- stiffened_gas_complete_safe.cpp ------------------===//
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
#include "materials/equations_of_state/stiffened_gas_complete_safe.h"
#include "utilities/helper_functions.h"
#include "utilities/mathematical_functions.h"
#include "utilities/string_operations.h"
#include <cmath>

/**
 * @brief Constructs a complete safe stiffened gas equation of state with eos
 * parameters given as input.
 * @param dimensional_eos_data Map containing all data for the equation of
 * state.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 *
 * @note During the constructing a check is done if the required parameter
 * exists. If not an error is thrown. Furthermore, dimensionalization of each
 * value is done.
 */
StiffenedGasCompleteSafe::StiffenedGasCompleteSafe(
    std::unordered_map<std::string, double> const &dimensional_eos_data,
    UnitHandler const &unit_handler)
    : gamma_(GetCheckedParameter<double>(dimensional_eos_data, "gamma",
                                         "StiffenedGasCompleteSafe")),
      energy_translation_factor_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_eos_data,
                                      "energyTranslationFactor",
                                      "StiffenedGasCompleteSafe"),
          {UnitType::Velocity, UnitType::Velocity}, {})),
      background_pressure_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_eos_data,
                                      "backgroundPressure",
                                      "StiffenedGasCompleteSafe"),
          UnitType::Pressure)),
      thermal_energy_factor_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_eos_data,
                                      "thermalEnergyFactor",
                                      "StiffenedGasCompleteSafe"),
          {UnitType::Velocity, UnitType::Velocity}, {UnitType::Temperature})),
      specific_gas_constant_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_eos_data,
                                      "specificGasConstant",
                                      "StiffenedGasCompleteSafe"),
          {UnitType::Velocity, UnitType::Velocity}, {UnitType::Temperature})) {
  /** Empty besides initializer list. */
}

/**
 * @brief Computes the pressure from inputs as -gamma * B + ( gamma - 1 ) * ( E
 * - 0.5 * ||v^2|| ) - rho * A.
 * @param mass The mass used for the computation.
 * @param momentum_x The momentum in x-direction used for the computation.
 * @param momentum_y The momentum in y-direction used for the computation.
 * @param momentum_z The momentum in z-direction used for the computation.
 * @param energy The energy used for the computation.
 * @return Pressure according to complete stiffened-gas equation of state.
 */
double StiffenedGasCompleteSafe::ComputePressure(double const mass,
                                                 double const momentum_x,
                                                 double const momentum_y,
                                                 double const momentum_z,
                                                 double const energy) const {
  double pressure =
      -gamma_ * background_pressure_ +
      (gamma_ - 1.0) *
          (energy -
           0.5 *
               DimensionAwareConsistencyManagedSum(momentum_x * momentum_x,
                                                   momentum_y * momentum_y,
                                                   momentum_z * momentum_z) /
               std::max(mass, epsilon_) -
           mass * energy_translation_factor_);
  return std::max(pressure, -background_pressure_ + epsilon_);
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
double StiffenedGasCompleteSafe::ComputeEnthalpy(double const mass,
                                                 double const momentum_x,
                                                 double const momentum_y,
                                                 double const momentum_z,
                                                 double const energy) const {
  return (energy +
          ComputePressure(mass, momentum_x, momentum_y, momentum_z, energy)) /
         std::max(mass, epsilon_);
}

/**
 * @brief Computes Energy according to stiffened gas complete equation, but
 * secured against division by zero.
 * @param density The density used for the computation.
 * @param velocity_x The velocity in x-direction used for the computation.
 * @param velocity_y The velocity in y-direction used for the computation.
 * @param velocity_z The velocity in z-direction used for the computation.
 * @param pressure The pressure used for the computation.
 * @return Energy according to given inputs.
 */
double StiffenedGasCompleteSafe::ComputeEnergy(double const density,
                                               double const velocity_x,
                                               double const velocity_y,
                                               double const velocity_z,
                                               double const pressure) const {
  return density * energy_translation_factor_ +
         (pressure + gamma_ * background_pressure_) / (gamma_ - 1.0) +
         (0.5 *
          DimensionAwareConsistencyManagedSum(velocity_x * velocity_x,
                                              velocity_y * velocity_y,
                                              velocity_z * velocity_z) *
          density);
}

/**
 * @brief Computes temperature for stiffened-gas EOS.
 * @param mass The mass used for the computation.
 * @param momentum_x The momentum in x-direction used for the computation.
 * @param momentum_y The momentum in y-direction used for the computation.
 * @param momentum_z The momentum in z-direction used for the computation.
 * @param energy The energy used for the computation.
 * @return Temperature according to complete stiffened-gas EOS.
 */
double StiffenedGasCompleteSafe::ComputeTemperature(double const mass,
                                                    double const momentum_x,
                                                    double const momentum_y,
                                                    double const momentum_z,
                                                    double const energy) const {
  return (ComputePressure(mass, momentum_x, momentum_y, momentum_z, energy) +
          background_pressure_) /
             (std::max(mass, epsilon_) * specific_gas_constant_) +
         ((gamma_ - 1.0) * thermal_energy_factor_ *
          std::pow(mass, gamma_ - 1.0)) /
             specific_gas_constant_;
}

/**
 * @brief Computes Gruneisen coefficient as ( gamma-1 ) for stiffened-gas
 * equation of state.
 * @return Gruneisen coefficient according to stiffened-gas equation of state.
 */
double StiffenedGasCompleteSafe::GetGruneisen() const { return (gamma_ - 1.0); }

/**
 * @brief Computes psi from inputs as ( p + gamma * B ) / rho.
 * @param pressure The pressure used for the computation.
 * @param one_density ( devision by zero is already avoided before ) .
 * @return Psi according to complete stiffened-gas equation of state.
 */
double StiffenedGasCompleteSafe::ComputePsi(double const pressure,
                                            double const one_density) const {
  return (pressure + gamma_ * background_pressure_) * one_density;
}

/**
 * @brief Returns Gamma.
 * @return Gamma.
 */
double StiffenedGasCompleteSafe::GetGamma() const { return gamma_; }

/**
 * @brief Returns B.
 * @return B.
 */
double StiffenedGasCompleteSafe::GetB() const { return background_pressure_; }

/**
 * @brief Computes Speed of Sound from inputs as sqrt( gamma * ( p + B ) / rho
 * ).
 * @param density The density used for the computation.
 * @param pressure The pressure used for the computation.
 * @return Speed of sound according to complete stiffened-gas equation of state.
 */
double
StiffenedGasCompleteSafe::ComputeSpeedOfSound(double const density,
                                              double const pressure) const {
  double speed_of_sound_squared =
      gamma_ * (pressure + background_pressure_) / std::max(density, epsilon_);
  return std::sqrt(std::max(speed_of_sound_squared, epsilon_));
}

/**
 * @brief Provides logging information of the equation of state.
 * @param indent Number of white spaces used at the beginning of each line for
 * the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string
StiffenedGasCompleteSafe::GetLogData(unsigned int const indent,
                                     UnitHandler const &unit_handler) const {
  // string that is returned
  std::string log_string;
  // Name of the equation of state
  log_string += StringOperations::Indent(indent) +
                "Type                       : Stiffened gas complete safe\n";
  // Parameters with small indentation
  log_string +=
      StringOperations::Indent(indent) + "Gruneisen coefficient      : " +
      StringOperations::ToScientificNotationString(GetGruneisen(), 9) + "\n";
  log_string += StringOperations::Indent(indent) +
                "Gamma                      : " +
                StringOperations::ToScientificNotationString(gamma_, 9) + "\n";
  log_string += StringOperations::Indent(indent) +
                "Energy translation factor  : " +
                StringOperations::ToScientificNotationString(
                    unit_handler.DimensionalizeValue(
                        energy_translation_factor_,
                        {UnitType::Velocity, UnitType::Velocity}, {}),
                    9) +
                "\n";
  log_string += StringOperations::Indent(indent) +
                "Stiffened pressure constant: " +
                StringOperations::ToScientificNotationString(
                    unit_handler.DimensionalizeValue(background_pressure_,
                                                     UnitType::Pressure),
                    9) +
                "\n";
  log_string +=
      StringOperations::Indent(indent) + "Thermal energy factor      : " +
      StringOperations::ToScientificNotationString(
          unit_handler.DimensionalizeValue(
              thermal_energy_factor_, {UnitType::Velocity, UnitType::Velocity},
              {UnitType::Temperature}),
          9) +
      "\n";
  log_string +=
      StringOperations::Indent(indent) + "Specific gas constant      : " +
      StringOperations::ToScientificNotationString(
          unit_handler.DimensionalizeValue(
              specific_gas_constant_, {UnitType::Velocity, UnitType::Velocity},
              {UnitType::Temperature}),
          9) +
      "\n";
  return log_string;
}
