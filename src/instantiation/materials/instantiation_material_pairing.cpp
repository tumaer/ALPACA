//===--------------- instantiation_material_pairing.cpp -------------------===//
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
#include "instantiation/materials/instantiation_material_pairing.h"

#include "materials/material_pairing_property_models/surface_tension_coefficient_models/constant_surface_tension_coefficient_model.h"
#include "materials/material_property_definitions.h"
#include "user_specifications/compile_time_constants.h"

namespace Instantiation {

/**
 * @brief Gives the model for the material pairing property surface tension
 * coefficient for the given input data.
 * @param model_name Name of the surface tension coefficient model.
 * @param model_data Data of the surface tension coefficient model.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @return pointer to the const base class of all material parameter models.
 */
std::unique_ptr<InterfaceParameterModel const>
InstantiateSurfaceTensionCoefficientModel(
    MaterialPropertyModelName const model_name,
    std::unordered_map<std::string, double> const &model_data,
    UnitHandler const &unit_handler) {
  // logger
  LogWriter &logger = LogWriter::Instance();

  // switch between different surface tension coefficient models
  switch (model_name) {
  case MaterialPropertyModelName::NotUsed: {
    return nullptr;
  }
  case MaterialPropertyModelName::SurfaceTensionCoefficientConstant: {
    // 1. Create, 2. Log, 3. Return model
    std::unique_ptr<ConstantSurfaceTensionCoefficientModel const> model(
        std::make_unique<ConstantSurfaceTensionCoefficientModel const>(
            model_data, unit_handler));
    logger.LogMessage(model->GetLogData(6, unit_handler));
    return model;
  }
  default: {
    throw std::logic_error(
        "This surface tension coefficient model is not implemented yet!");
  }
  }
}

/**
 * @brief Instantiates a complete material pairing with the given input reader.
 * @param material_indices Indices of the materials that form the pairing, which
 * should be read from the input file.
 * @param material_reader Reader that provides access to the material data of
 * the input file.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @return Full initialized material pairing.
 *
 * @note during the initialization checks are done to read only variables that
 * are required with the given compile time settings.
 */
MaterialPairing
InstantiateMaterialPairing(std::vector<unsigned int> const &material_indices,
                           MaterialReader const &material_reader,
                           UnitHandler const &unit_handler) {

  // Read all data from input file
  // Declare all variables that can be filled (use default values)
  double surface_tension_coefficient_fixed_value = -1.0;
  MaterialPropertyModelName surface_tension_coefficient_model_name =
      MaterialPropertyModelName::NotUsed;
  std::unordered_map<std::string, double>
      surface_tension_coefficient_model_data = {};
  std::unique_ptr<InterfaceParameterModel const>
      surface_tension_coefficient_model = nullptr;

  // Surface tension coefficient only if needed
  if constexpr (CC::CapillaryForcesActive()) {
    // read model data if required
    if constexpr (CC::SurfaceTensionCoefficientModelActive()) {
      // read the model data
      surface_tension_coefficient_model_name = material_reader.ReadModelName(
          material_indices, MaterialProperty::SurfaceTensionCoefficient);
      surface_tension_coefficient_model_data = material_reader.ReadModelData(
          material_indices, MaterialProperty::SurfaceTensionCoefficient);
    } else { // otherwise read the fixed value
      surface_tension_coefficient_fixed_value = material_reader.ReadFixedValue(
          material_indices, MaterialProperty::SurfaceTensionCoefficient);
    }
  }

  // Create final data (eos + modles) and log information
  // logger
  LogWriter &logger = LogWriter::Instance();
  logger.LogMessage(" ");
  logger.LogMessage("Material pairing " + std::to_string(material_indices[0]) +
                    " <-> " + std::to_string(material_indices[1]) + ":");

  // Material pairing properties
  // Log properties only if positive values are specified (they exist)
  std::string tmp_string;
  // Bulk viscosity
  tmp_string = StringOperations::Indent(2) + "Surface tension coefficient: ";
  if constexpr (CC::CapillaryForcesActive()) {
    if (surface_tension_coefficient_model_name !=
        MaterialPropertyModelName::NotUsed) {
      // create model (logs data itself)
      logger.LogMessage(tmp_string);
      surface_tension_coefficient_model =
          InstantiateSurfaceTensionCoefficientModel(
              surface_tension_coefficient_model_name,
              surface_tension_coefficient_model_data, unit_handler);
    } else {
      logger.LogMessage(tmp_string +
                        StringOperations::ToScientificNotationString(
                            surface_tension_coefficient_fixed_value, 9));
    }
  } else {
    logger.LogMessage(tmp_string + "Not Used");
  }

  // return the fully initialied material pairing
  return MaterialPairing(surface_tension_coefficient_fixed_value,
                         std::move(surface_tension_coefficient_model),
                         unit_handler);
}

} // namespace Instantiation
