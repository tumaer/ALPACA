//===-------------- constant_shear_viscosity_model.cpp --------------------===//
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
#include "materials/material_property_models/shear_viscosity_models/constant_shear_viscosity_model.h"

#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"

/**
 * @brief Constructor to instantiate the constant shear viscosity model.
 * @param dimensional_parameter_map Map, where all parameters are stored.
 * @param unit_handler Instance to provide dimensionalization of values.
 *
 * @note Runtime a check is done that all required parameters are present.
 * Furthermore, pre-calculations are done for simpler access during compute
 * call.
 */
ConstantShearViscosityModel::ConstantShearViscosityModel(
    std::unordered_map<std::string, double> const &dimensional_parameter_map,
    UnitHandler const &unit_handler)
    : // Start initializer list
      ConstantMaterialParameterModel<ConstantShearViscosityModel>(),
      mu_constant_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_parameter_map, "muConstant",
                                      "ConstantShearViscosityModel"),
          UnitType::Viscosity)) {
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
ConstantShearViscosityModel::GetLogData(unsigned int const indent,
                                        UnitHandler const &unit_handler) const {
  // string that is returned
  std::string log_string;
  // Add data
  log_string += StringOperations::Indent(indent) + "Model type: Constant \n";
  log_string +=
      StringOperations::Indent(indent) + "mu_const  : " +
      StringOperations::ToScientificNotationString(
          unit_handler.DimensionalizeValue(mu_constant_, UnitType::Viscosity),
          9) +
      "\n";

  return log_string;
}

/**
 * @brief Provides a constant shear viscosity from user input
 * @return The constant shear viscosity
 */
double ConstantShearViscosityModel::ComputeParameter() const {
  return mu_constant_;
}
