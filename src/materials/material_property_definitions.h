//===--------------- material_property_definitions.h ----------------------===//
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
#ifndef MATERIAL_PROPERTY_DEFINITIONS_H
#define MATERIAL_PROPERTY_DEFINITIONS_H

#include <stdexcept>
#include <string>

#include "utilities/string_operations.h"

/**
 * @brief The MaterialPropertyName defines the identifiers to all different pure
 * material and material pairing properties that exist. It is used to provide a
 * defined mapping for the properties between input, transfer and output.
 * @note The identifier is also used for the mapping of specific property
 * parameter models.
 */
enum class MaterialProperty {
  ShearViscosity,
  BulkViscosity,
  ThermalConductivity,
  SpecificHeatCapacity,
  SurfaceTensionCoefficient
};

/**
 * @brief Converts the material property to a proper string that can be used for
 * input and output of the data.
 * @param property The material property identifier.
 * @return first_capitalized Flag to indicate if the first letter should be
 * capitalized or not.
 *
 * @note Do not change the string obtained for the given property, since this is
 * also used for mapping the property to a corresponding parameter model. In
 * case change both.
 */
inline std::string
MaterialPropertyToString(MaterialProperty const property,
                         bool const first_capitalized = false) {

  switch (property) {
  case MaterialProperty::SpecificHeatCapacity: {
    return first_capitalized ? "SpecificHeatCapacity" : "specificHeatCapacity";
  }
  case MaterialProperty::ThermalConductivity: {
    return first_capitalized ? "ThermalConductivity" : "thermalConductivity";
  }
  case MaterialProperty::ShearViscosity: {
    return first_capitalized ? "ShearViscosity" : "shearViscosity";
  }
  case MaterialProperty::BulkViscosity: {
    return first_capitalized ? "BulkViscosity" : "bulkViscosity";
  }
  case MaterialProperty::SurfaceTensionCoefficient: {
    return first_capitalized ? "SurfaceTensionCoefficient"
                             : "surfaceTensionCoefficient";
  }
  default: {
    // if nothing matches throw error
    throw std::logic_error("Material property not known!");
  }
  }
}

/**
 * @brief The MaterialPropertyModelName class defines all parameter models
 * designed for computing material properties on other given material quantities
 *        such as prime states.
 *
 * @note The name MUST follow the syntax material property + model (exception is
 * NotUsed).
 */
enum class MaterialPropertyModelName {
  // Default name if model is not used
  NotUsed,
  // viscosity models
  // constant model
  ShearViscosityConstant,
  // Shear rate models
  ShearViscosityPowerLaw,
  ShearViscosityCross,
  ShearViscosityCarreauYasuda,

  // Temperature models
  ShearViscositySutherland,

  // Thermal conductivity models
  // constant model
  ThermalConductivityConstant,

  // Surface tension coefficient models
  // constant model
  SurfaceTensionCoefficientConstant
};

/**
 * @brief Converts a string to its corresponding MaterialPropertyModelName.
 * @param property Material property for which the model is used (e.g.
 * ShearViscosity, ThermalConductivity).
 * @param model_name Name of the model that is desired (e.g., Constant,
 * Carreau).
 * @return Name of the material property model.
 */
inline MaterialPropertyModelName
StringToMaterialPropertyModel(MaterialProperty const property,
                              std::string const &model_name) {
  // transform string to upper case without spaces
  std::string const name_upper_case(StringOperations::ToUpperCaseWithoutSpaces(
      MaterialPropertyToString(property) + model_name));
  // Viscosity models
  // constant
  if (name_upper_case == "SHEARVISCOSITYCONSTANT") {
    return MaterialPropertyModelName::ShearViscosityConstant;
  }
  // shear rate models
  else if (name_upper_case == "SHEARVISCOSITYPOWERLAW") {
    return MaterialPropertyModelName::ShearViscosityPowerLaw;
  } else if (name_upper_case == "SHEARVISCOSITYCROSS") {
    return MaterialPropertyModelName::ShearViscosityCross;
  } else if (name_upper_case == "SHEARVISCOSITYCARREAUYASUDA") {
    return MaterialPropertyModelName::ShearViscosityCarreauYasuda;
  }
  // temperature models
  else if (name_upper_case == "SHEARVISCOSITYSUTHERLAND") {
    return MaterialPropertyModelName::ShearViscositySutherland;
  }

  // conductivity models
  // constant
  else if (name_upper_case == "THERMALCONDUCTIVITYCONSTANT") {
    return MaterialPropertyModelName::ThermalConductivityConstant;
  }

  // surface tension coefficient models
  // constant
  else if (name_upper_case == "SURFACETENSIONCOEFFICIENTCONSTANT") {
    return MaterialPropertyModelName::SurfaceTensionCoefficientConstant;
  }

  // default behavior if nothing is known
  else {
    throw std::logic_error("Material property model is not known!");
  }
}

/**
 * @brief Converts a MaterialPropertyModelName into a corresponding string
 * @param name The identifier of the material property model
 * @return String to be used
 */
inline std::string
MaterialPropertyModelToString(MaterialPropertyModelName const name) {

  switch (name) {
  // All constant models
  case MaterialPropertyModelName::ShearViscosityConstant:
  case MaterialPropertyModelName::ThermalConductivityConstant:
  case MaterialPropertyModelName::SurfaceTensionCoefficientConstant: {
    return "Constant";
  }
  // All specific models
  case MaterialPropertyModelName::ShearViscosityPowerLaw: {
    return "PowerLaw";
  }
  case MaterialPropertyModelName::ShearViscosityCross: {
    return "Cross";
  }
  case MaterialPropertyModelName::ShearViscosityCarreauYasuda: {
    return "CarreauYasuda";
  }
  case MaterialPropertyModelName::ShearViscositySutherland: {
    return "Sutherland";
  }
  // Default model not used
  default: {
    return "Not Used";
  }
  }
}

#endif // MATERIAL_PROPERTY_DEFINITIONS_H
