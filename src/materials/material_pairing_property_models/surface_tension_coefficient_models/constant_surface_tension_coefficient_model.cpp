//===-------- constant_surface_tension_coefficient_model.cpp --------------===//
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
#include "materials/material_pairing_property_models/surface_tension_coefficient_models/constant_surface_tension_coefficient_model.h"
#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"

/**
 * @brief Constructor to instantiate the constant thermal conductivity model
 * @param parameter_map Map, where all parameters are stored
 * @param unit_handler Instance to provide dimensionalization of values
 *
 * @note Runtime a check is done that all required parameters are present.
 * Furthermore, pre-calculations are done to simply access in compute call
 */
ConstantSurfaceTensionCoefficientModel::ConstantSurfaceTensionCoefficientModel(
    std::unordered_map<std::string, double> const &parameter_map,
    UnitHandler const &unit_handler)
    : // Start initializer list
      ConstantInterfaceParameterModel<ConstantSurfaceTensionCoefficientModel>(),
      sigma_constant_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(parameter_map, "sigmaConstant",
                                      "ConstantSurfaceTensionCoefficientModel"),
          UnitType::SurfaceTensionCoefficient)) {
  /** Empty besides initializer list and friend class constructor call  */
}

/**
 * @brief Provides logging information of the given model.
 * @param log_string String that is filled with information.
 */
std::string ConstantSurfaceTensionCoefficientModel::GetLogData(
    unsigned int const indent, UnitHandler const &unit_handler) const {
  // string that is returned
  std::string log_string;
  // Add data
  log_string += StringOperations::Indent(indent) + "Model type: Constant \n";
  log_string += StringOperations::Indent(indent) + "sigma0    : " +
                StringOperations::ToScientificNotationString(
                    unit_handler.DimensionalizeValue(
                        sigma_constant_, UnitType::SurfaceTensionCoefficient),
                    9) +
                "\n";

  return log_string;
}

/**
 * @brief Provides a constant surface tension coefficient from user input.
 * @return The constant surface tension coefficient.
 */
double ConstantSurfaceTensionCoefficientModel::ComputeParameter() const {
  return sigma_constant_;
}
