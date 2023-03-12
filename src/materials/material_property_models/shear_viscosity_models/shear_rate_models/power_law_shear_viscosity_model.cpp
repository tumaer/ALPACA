//===--------------- power_law_shear_viscosity_model.cpp ------------------===//
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
#include "materials/material_property_models/shear_viscosity_models/shear_rate_models/power_law_shear_viscosity_model.h"

#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"

/**
 * @brief Constructor to instantiate the Power Law shear viscosity model. The
 * model is implemented based on \cite{Deville2012}. K         = consistency
 * factor [unit depends on n] (for n != 1), Newtonian viscosity for n = 1 mu = K
 * * (gamma_dot)^(n-1)       gamma_dot   = shear rate n      = power law index
 * (n=1: Newtonian, n<1: shear thinning, n>1: shear thickening/dilating)
 *
 * @param dimensional_parameter_map Map, where all parameters are stored.
 * @param unit_handler Instance to provide dimensionalization of values.
 *
 * @note Runtime a check is done that all required parameters are present.
 * Furthermore, pre-calculations are done for simpler access during compute call
 */
PowerLawShearViscosityModel::PowerLawShearViscosityModel(
    std::unordered_map<std::string, double> const &dimensional_parameter_map,
    UnitHandler const &unit_handler)
    : // Start initializer list
      ShearRateMaterialParameterModel<PowerLawShearViscosityModel>(),
      consistency_factor_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_parameter_map,
                                      "consistencyFactor",
                                      "PowerLawShearViscosityModel"),
          UnitType::Viscosity)),
      power_law_exponent_(GetCheckedParameter<double>(
          dimensional_parameter_map, "powerLawExponent",
          "PowerLawShearViscosityModel")),
      exponent_(power_law_exponent_ - 1.0) {
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
PowerLawShearViscosityModel::GetLogData(unsigned int const indent,
                                        UnitHandler const &unit_handler) const {
  // string that is returned
  std::string log_string;
  // Add data
  log_string += StringOperations::Indent(indent) + "Model type: Power-law \n";
  log_string += StringOperations::Indent(indent) + "factor    : " +
                StringOperations::ToScientificNotationString(
                    unit_handler.DimensionalizeValue(consistency_factor_,
                                                     UnitType::Viscosity),
                    9) +
                "\n";
  log_string +=
      StringOperations::Indent(indent) + "n         : " +
      StringOperations::ToScientificNotationString(power_law_exponent_, 9) +
      "\n";

  return log_string;
}

/**
 * @brief Actual implementation of the shear rate model.
 * @param shear_rate The shear rate for which the viscosity is computed.
 * @return The calculated shear viscosity.
 */
double
PowerLawShearViscosityModel::ComputeParameter(double const shear_rate) const {
  return consistency_factor_ * std::pow(shear_rate, exponent_);
}
