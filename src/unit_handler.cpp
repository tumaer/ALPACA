//===----------------------- unit_handler.cpp -----------------------------===//
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
#include "unit_handler.h"

#include <stdexcept>
#include <string>

#include "utilities/string_operations.h"

/**
 * @brief Default constructor.
 * @param parser A parser for the user input file.
 */
UnitHandler::UnitHandler(double const density_reference,
                         double const velocity_reference,
                         double const length_reference,
                         double const temperature_reference)
    : density_reference_(density_reference),
      velocity_reference_(velocity_reference),
      length_reference_(length_reference),
      temperature_reference_(temperature_reference),
      time_reference_(length_reference_ / velocity_reference_),
      mass_reference_(density_reference_ * length_reference_ *
                      length_reference_ * length_reference_),
      momentum_reference_(density_reference_ * velocity_reference_),
      pressure_reference_(density_reference_ * velocity_reference_ *
                          velocity_reference_),
      energy_reference_(density_reference_ * velocity_reference_ *
                        velocity_reference_),
      viscosity_reference_(density_reference_ * velocity_reference_ *
                           length_reference_),
      thermal_conductivity_reference_(
          density_reference_ * velocity_reference_ * velocity_reference_ *
          velocity_reference_ * length_reference_ / temperature_reference_),
      surface_tension_coefficient_reference_(
          density_reference_ * velocity_reference_ * velocity_reference_ *
          length_reference_) {
  /** Empty besides initializer list */
}

/**
 * @brief Translates a unit based value into a unitless one, using the specific
 * function calls of the given unit type.
 * @param value value which should be made dimensionless.
 * @param unit_type Unit of the value..
 * @return Unitless value.
 */
double UnitHandler::NonDimensionalizeValue(double const value,
                                           UnitType const unit_type) const {
  switch (unit_type) {
  case UnitType::Unitless: {
    return value;
  }
  case UnitType::Length: {
    return value / length_reference_;
  }
  case UnitType::Time: {
    return value / time_reference_;
  }
  case UnitType::Mass: {
    return value / mass_reference_;
  }
  case UnitType::Density: {
    return value / density_reference_;
  }
  case UnitType::Velocity: {
    return value / velocity_reference_;
  }
  case UnitType::Momentum: {
    return value / momentum_reference_;
  }
  case UnitType::Pressure: {
    return value / pressure_reference_;
  }
  case UnitType::Temperature: {
    return value / temperature_reference_;
  }
  case UnitType::Energy: {
    return value / energy_reference_;
  }
  case UnitType::Viscosity: {
    return value / viscosity_reference_;
  }
  case UnitType::ThermalConductivity: {
    return value / thermal_conductivity_reference_;
  }
  case UnitType::SurfaceTensionCoefficient: {
    return value / surface_tension_coefficient_reference_;
  }
  default: {
    throw std::logic_error("UnitType " +
                           std::to_string(static_cast<int>(unit_type)) +
                           " is unknown and cannot be non-dimensionalized");
  } // suppresses Compiler Warning "control reaches end of non-void function
    // [-Wreturn-type]"
  }
}

/**
 * @brief Translates a unitless value into a dimensional one, using the specific
 * function calls of the given unit type.
 * @param value value which should be dimensionalized.
 * @param unit_type Unit of the value.
 * @return Unitless value.
 */
double UnitHandler::DimensionalizeValue(double const value,
                                        UnitType const unit_type) const {
  switch (unit_type) {
  case UnitType::Unitless: {
    return value;
  }
  case UnitType::Length: {
    return value * length_reference_;
  }
  case UnitType::Time: {
    return value * time_reference_;
  }
  case UnitType::Mass: {
    return value * mass_reference_;
  }
  case UnitType::Density: {
    return value * density_reference_;
  }
  case UnitType::Velocity: {
    return value * velocity_reference_;
  }
  case UnitType::Momentum: {
    return value * momentum_reference_;
  }
  case UnitType::Pressure: {
    return value * pressure_reference_;
  }
  case UnitType::Temperature: {
    return value * temperature_reference_;
  }
  case UnitType::Energy: {
    return value * energy_reference_;
  }
  case UnitType::Viscosity: {
    return value * viscosity_reference_;
  }
  case UnitType::ThermalConductivity: {
    return value * thermal_conductivity_reference_;
  }
  case UnitType::SurfaceTensionCoefficient: {
    return value * surface_tension_coefficient_reference_;
  }
  default: {
    throw std::logic_error("UnitType " +
                           std::to_string(static_cast<int>(unit_type)) +
                           " is unknown and cannot be dimensionalized");
  } // suppresses Compiler Warning "control reaches end of non-void function
    // [-Wreturn-type]"
  }
}

/**
 * @brief Non-dimensionalizes a given value with a set of nominator and
 * denominator units.
 * @param value Dimensional value that should be non-dimensionalized.
 * @param nominator_units All units used in the nominator.
 * @param denominator_units All units used in the denominator.
 * @return Dimensionless value.
 */
double UnitHandler::NonDimensionalizeValue(
    double const value, std::vector<UnitType> const &nominator_units,
    std::vector<UnitType> const &denominator_units) const {
  // declare nominator and denominator values
  double nominator = 1.0;
  double denominator = 1.0;
  // Loop through all nominator and get nominator value (use dimensionalize
  // function since the multiplication is required)
  for (UnitType const &unit_type : nominator_units) {
    nominator *= DimensionalizeValue(1.0, unit_type);
  }
  // Loop through all denominator units and get denominator value (use
  // dimensionalize function since the multiplication is required)
  for (UnitType const &unit_type : denominator_units) {
    denominator *= DimensionalizeValue(1.0, unit_type);
  }
  // reverse order required for non-dimensionalization
  return value * denominator / nominator;
}

/**
 * @brief Dimensionalizes a given value with a set of nominator and denominator
 * units.
 * @param value Dimensional value that should be dimensionalized.
 * @param nominator_units All units used in the nominator.
 * @param denominator_units All units used in the denominator.
 * @return Dimensional value.
 */
double UnitHandler::DimensionalizeValue(
    double const value, std::vector<UnitType> const &nominator_units,
    std::vector<UnitType> const &denominator_units) const {
  // declare nominator and denominator values
  double nominator = 1.0;
  double denominator = 1.0;
  // Loop through all nominator and get nominator value (use dimensionalize
  // function since the multiplication is required)
  for (UnitType const &unit_type : nominator_units) {
    nominator *= DimensionalizeValue(1.0, unit_type);
  }
  // Loop through all denominator units and get denominator value (use
  // dimensionalize function since the multiplication is required)
  for (UnitType const &unit_type : denominator_units) {
    denominator *= DimensionalizeValue(1.0, unit_type);
  }
  return value * nominator / denominator;
}
