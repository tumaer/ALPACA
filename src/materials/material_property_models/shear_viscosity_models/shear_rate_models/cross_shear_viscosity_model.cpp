//===--------------- cross_shear_viscosity_model.cpp ----------------------===//
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
#include "materials/material_property_models/shear_viscosity_models/shear_rate_models/cross_shear_viscosity_model.h"

#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"

/**
 * @brief Constructor to instantiate the Cross shear viscosity model. The model
 * is implemented based on \cite{Cross1965}. mu0  - mu_inf          mu0/mu-inf
 * = viscosity at zero/infinite shear rates. mu = mu_inf +
 * -------------------------        k       = the shear rate where the viscosity
 * is half between mu0 and muInf. 1 + ( gamma_dot / k )^n     gamma_dot = shear
 * rate. n       = power-law index (no recovering of Newtonian behavior
 * possible, but n=0 gives constant shear viscosity at half between mu0 and
 * mu_inf).
 *
 * @param dimensional_parameter_map Map, where all parameters are stored.
 * @param unit_handler Instance to provide dimensionalization of values.
 *
 * @note Runtime a check is done that all required parameters are present.
 * Furthermore, pre-calculations are done for simpler access during compute
 * call.
 */
CrossShearViscosityModel::CrossShearViscosityModel(
    std::unordered_map<std::string, double> const &dimensional_parameter_map,
    UnitHandler const &unit_handler)
    : // Start initializer list
      ShearRateMaterialParameterModel<CrossShearViscosityModel>(),
      mu_infinite_shear_rates_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_parameter_map,
                                      "muInfiniteShearRates",
                                      "CrossShearViscosityModel"),
          UnitType::Viscosity)),
      mu_zero_shear_rates_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_parameter_map,
                                      "muZeroShearRates",
                                      "CrossShearViscosityModel"),
          UnitType::Viscosity)),
      power_law_exponent_(GetCheckedParameter<double>(
          dimensional_parameter_map, "powerLawExponent",
          "CrossShearViscosityModel")),
      shear_rate_mu_half_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_parameter_map,
                                      "shearRateHalfMu",
                                      "CrossShearViscosityModel"),
          {}, {UnitType::Time})),
      mu_zero_minus_infinite_(mu_zero_shear_rates_ - mu_infinite_shear_rates_),
      one_shear_rate_mu_half_(1.0 / shear_rate_mu_half_) {
  /** Empty besides initializer list and friend class constructor call  */
}

/**
 * @brief Provides logging information of the given model.
 * @param indent Number of white spaces used at the beginning of each line for
 * the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string
CrossShearViscosityModel::GetLogData(unsigned int const indent,
                                     UnitHandler const &unit_handler) const {
  // string that is returned
  std::string log_string;
  // Add data
  log_string += StringOperations::Indent(indent) + "Model type: Cross \n";
  log_string += StringOperations::Indent(indent) + "mu_inf    : " +
                StringOperations::ToScientificNotationString(
                    unit_handler.DimensionalizeValue(mu_infinite_shear_rates_,
                                                     UnitType::Viscosity),
                    9) +
                "\n";
  log_string += StringOperations::Indent(indent) + "mu0       : " +
                StringOperations::ToScientificNotationString(
                    unit_handler.DimensionalizeValue(mu_zero_shear_rates_,
                                                     UnitType::Viscosity),
                    9) +
                "\n";
  log_string +=
      StringOperations::Indent(indent) + "n         : " +
      StringOperations::ToScientificNotationString(power_law_exponent_, 9) +
      "\n";
  log_string +=
      StringOperations::Indent(indent) + "k         : " +
      StringOperations::ToScientificNotationString(shear_rate_mu_half_, 9) +
      "\n";

  return log_string;
}

/**
 * @brief Actual implementation of the shear rate model.
 * @param shear_rate The shear rate for which the viscosity is computed.
 * @return The calculated shear viscosity.
 */
double
CrossShearViscosityModel::ComputeParameter(double const shear_rate) const {
  double const denominator =
      1.0 + std::pow(shear_rate * one_shear_rate_mu_half_, power_law_exponent_);
  return mu_infinite_shear_rates_ + mu_zero_minus_infinite_ / denominator;
}
