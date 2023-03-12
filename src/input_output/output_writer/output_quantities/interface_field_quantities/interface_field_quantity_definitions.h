//===------------- interface_field_quantity_definitions.h -----------------===//
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
#ifndef INTERFACE_FIELD_QUANTITIES_H
#define INTERFACE_FIELD_QUANTITIES_H

#include "block_definitions/field_details.h"
#include "block_definitions/field_interface_definitions.h"
#include "enums/interface_tag_definition.h"
#include <iostream>
#include <vector>

/**
 * @brief The InterfaceFieldName enum gives all possible choices for which an
 * output can be written (interface descriptions, interface states and interface
 * paramteres).
 *
 * @note: In order to append a InterfaceFieldQuantity to the ouput follow the
 * steps:
 *        1. Add The field quantity to InterfaceFieldQuantityName in the
 * appropriate group (arbitrary name)
 *        2. Add for this name an entry with default values to the function
 * DataOfInterfaceFieldQuantity below ( for scalar quantities follow e.g.
 * levelset, for vectorial quantities follow e.g. velocity in
 * material_field_quantity_definitions.h)
 *        3. Add an output constants expression in
 * user_specifications/output_constants.h
 *        4. Add the constructor call in the "instantiation_output_writer.cpp"
 * function.
 */
enum class InterfaceFieldQuantityName {
  // interface descriptions
  Levelset,
  VolumeFraction,
  // interface states
  InterfaceVelocity,
  PressurePositive,
  PressureNegative,
  // interface parameters
  SurfaceTensionCoefficient
};

/**
 * @brief Returns the suffix for a given interface buffer type.
 * @param buffer_type interface buffer type identifier.
 * @return suffix for the given identifier.
 */
inline std::string SuffixOfInterfaceDescriptionBufferType(
    InterfaceDescriptionBufferType const buffer_type) {
  switch (buffer_type) {
  case InterfaceDescriptionBufferType::Base: {
    return "";
  }
  case InterfaceDescriptionBufferType::RightHandSide: {
    return "_rhs";
  }
  case InterfaceDescriptionBufferType::Reinitialized: {
    return "_reinitialized";
  }
  default: { // last possibility InterfaceDescriptionBufferType::Initial :
    return "_initial";
  }
  }
}

/**
 * @brief Defines a struct storing all necessary information to provide the
 * correct information for each InterfaceFieldQuantity.
 */
struct InterfaceFieldQuantityData {
  // Interface field type for which the output is written
  InterfaceFieldType field_type_;
  // indices of the field type which should be included in the output (more than
  // one element generates vectorial output)
  std::vector<unsigned int> field_indices_;
  // flag wether interface tags should be used for the output
  bool use_interface_tags_for_default_;
  // default value for the standard and debug output
  double default_value_;

  /**
   * @brief Gives the dimensions of the quantity.
   * @return Dimensions of the quantity.
   */
  std::array<unsigned int, 2> GetDimensions() const {
    return {static_cast<unsigned int>(field_indices_.size()), 1};
  }

  /**
   * @brief Verifies the quantity data on validity to eliminate non-active
   * fields.
   * @return True if any of the fields are active, otherwise false.
   * @note Manipulates the member data in case some fields are not active.
   */
  bool VerifyQuantity() {
    // Local field type for lambda capture
    InterfaceFieldType const field_type = field_type_;
    // Store current size for warning handling
    size_t current_size = field_indices_.size();
    // Erase all elements that are not active
    field_indices_.erase(
        std::remove_if(field_indices_.begin(), field_indices_.end(),
                       [&field_type](double const &index) {
                         return !IF::IsFieldActive(field_type, index);
                       }),
        field_indices_.end());
    // Throw warning if some of the fields have been rejected
    if (current_size != field_indices_.size()) {
      std::cerr << "Output Warning! Some of the desired interface fields are "
                   "not active! Continue with reduced set!"
                << std::endl;
    }
    // If no field index is left anymore return false otherwise true
    return !field_indices_.empty();
  }
};

/***********************************************************************************************************************/
/*                              Add here new interface field quantities */
/***********************************************************************************************************************/
/**
 * @brief Gives the appropriate dat for the interface field quantity.
 * @param output_quantity The Interface field output quantity identifier.
 * @return Interface Field Type of given identifier.
 */
inline InterfaceFieldQuantityData
DataOfInterfaceFieldQuantity(InterfaceFieldQuantityName const output_quantity) {
  // General declaration of the struct
  InterfaceFieldQuantityData quantity_data;
  // Differ between all quantities
  switch (output_quantity) {
  // interface descriptions
  case InterfaceFieldQuantityName::Levelset: {
    quantity_data.field_type_ = InterfaceFieldType::Description;
    quantity_data.field_indices_ = {IDTI(InterfaceDescription::Levelset)};
    quantity_data.use_interface_tags_for_default_ = true;
    quantity_data.default_value_ = CC::LSCOF();
    return quantity_data;
  }
  case InterfaceFieldQuantityName::VolumeFraction: {
    quantity_data.field_type_ = InterfaceFieldType::Description;
    quantity_data.field_indices_ = {IDTI(InterfaceDescription::VolumeFraction)};
    quantity_data.use_interface_tags_for_default_ = true;
    quantity_data.default_value_ = 1.0 / ITTI(IT::BulkPhase);
    return quantity_data;
  }
  // interface states
  case InterfaceFieldQuantityName::InterfaceVelocity: {
    quantity_data.field_type_ = InterfaceFieldType::States;
    quantity_data.field_indices_ = {ISTI(InterfaceState::Velocity)};
    quantity_data.use_interface_tags_for_default_ = false;
    quantity_data.default_value_ = 0.0;
    return quantity_data;
  }
  case InterfaceFieldQuantityName::PressurePositive: {
    quantity_data.field_type_ = InterfaceFieldType::States;
    quantity_data.field_indices_ = {ISTI(InterfaceState::PressurePositive)};
    quantity_data.use_interface_tags_for_default_ = false;
    quantity_data.default_value_ = 0.0;
    return quantity_data;
  }
  case InterfaceFieldQuantityName::PressureNegative: {
    quantity_data.field_type_ = InterfaceFieldType::States;
    quantity_data.field_indices_ = {ISTI(InterfaceState::PressureNegative)};
    quantity_data.use_interface_tags_for_default_ = false;
    quantity_data.default_value_ = 0.0;
    return quantity_data;
  }
  // interface parameters
  case InterfaceFieldQuantityName::SurfaceTensionCoefficient: {
    quantity_data.field_type_ = InterfaceFieldType::Parameters;
    quantity_data.field_indices_ = {
        IPTI(InterfaceParameter::SurfaceTensionCoefficient)};
    quantity_data.use_interface_tags_for_default_ = false;
    quantity_data.default_value_ = 0.0;
    return quantity_data;
  }
  // default if nothing matches
  default: {
    throw std::logic_error("Interface field quantity not known!");
  }
  }
}

#endif // INTERFACE_FIELD_QUANTITIES_H
