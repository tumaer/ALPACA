//===------------- carreau_yasuda_shear_viscosity_model.cpp ---------=-----===//
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
#include "materials/material_property_models/shear_viscosity_models/shear_rate_models/carreau_yasuda_shear_viscosity_model.h"

#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"

/**
 * @brief Constructor to instantiate the Carreau-Yasuda shear viscosity model.
 * The model is implemented based on \cite{Bird1987}. mu0/mu-inf = viscosity at
 * zero/infinite shear rates. mu = mu_inf + (mu0  - mu_inf) * (1 + (gamma_dot /
 * k)^a)^((n-1)/a)          k      = the shear rate where the viscosity is
 * constant (depending on n and a). gamma_dot = shear rate. n      = power-law
 * exponent (slope in power law region). a      = transition factor between
 * zero-shear-rate and power-law region. A transition factor of 2 recovers the
 * original Carreau model.
 *
 * @param dimensional_parameter_map Map, where all parameters are stored.
 * @param unit_handler Instance to provide dimensionalization of values.
 *
 * @note Runtime a check is done that all required parameters are present.
 * Furthermore, pre-calculations are done for simpler access during compute
 * call.
 */
CarreauYasudaShearViscosityModel::CarreauYasudaShearViscosityModel(
    std::unordered_map<std::string, double> const &dimensional_parameter_map,
    UnitHandler const &unit_handler)
    : // Start initializer list
      ShearRateMaterialParameterModel<CarreauYasudaShearViscosityModel>(),
      mu_infinite_shear_rates_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_parameter_map,
                                      "muInfiniteShearRates",
                                      "CarreauYasudaShearViscosityModel"),
          UnitType::Viscosity)),
      mu_zero_shear_rates_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_parameter_map,
                                      "muZeroShearRates",
                                      "CarreauYasudaShearViscosityModel"),
          UnitType::Viscosity)),
      transition_factor_(GetCheckedParameter<double>(
          dimensional_parameter_map, "transitionFactor",
          "CarreauYasudaShearViscosityModel")),
      power_law_exponent_(GetCheckedParameter<double>(
          dimensional_parameter_map, "powerLawExponent",
          "CarreauYasudaShearViscosityModel")),
      shear_rate_constant_mu_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_parameter_map,
                                      "shearRateConstantMu",
                                      "CarreauYasudaShearViscosityModel"),
          {}, {UnitType::Time})),
      mu_zero_minus_infinite_(mu_zero_shear_rates_ - mu_infinite_shear_rates_),
      exponent_((power_law_exponent_ - 1.0) / transition_factor_),
      one_shear_rate_constant_mu_(1.0 / shear_rate_constant_mu_) {
  /** Empty besides initializer list and friend class constructor call  */
}

/**
 * @brief Provides logging information of the given model.
 * @param indent Number of white spaces used at the beginning of each line for
 * the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string CarreauYasudaShearViscosityModel::GetLogData(
    unsigned int const indent, UnitHandler const &unit_handler) const {
  // string that is returned
  std::string log_string;
  // Add data
  log_string +=
      StringOperations::Indent(indent) + "Model type: Carreau-Yasuda \n";
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
      StringOperations::Indent(indent) + "a         : " +
      StringOperations::ToScientificNotationString(transition_factor_, 9) +
      "\n";
  log_string +=
      StringOperations::Indent(indent) + "constant  : " +
      StringOperations::ToScientificNotationString(shear_rate_constant_mu_, 9) +
      "\n";

  return log_string;
}

/**
 * @brief Computes the shear viscosity based on the model parameter and given
 * shear rate.
 * @param shear_rate The shear rate for which the viscosity is computed.
 * @return The calculated shear viscosity.
 */
double CarreauYasudaShearViscosityModel::ComputeParameter(
    double const shear_rate) const {

  double const base = 1.0 + std::pow(shear_rate * one_shear_rate_constant_mu_,
                                     transition_factor_);
  return mu_infinite_shear_rates_ +
         mu_zero_minus_infinite_ * std::pow(base, exponent_);
}
