//===------------------ interface_field_quantity.cpp ----------------------===//
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
#include "input_output/output_writer/output_quantities/interface_field_quantities/interface_field_quantity.h"
#include "input_output/utilities/xdmf_utilities.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief constructor to create interface field output.
 * @param unit_handler Unit handler class for dimensionalization.
 * @param material_manager Material manager for accessing material data.
 * @param quantity_name The name of the interface field quantity used for the
 * output.
 * @param output_flags Flags for which output type an output is written (0:
 * standard, 1: interface, 2:debug).
 * @param quantity_data Data struct that contains all relevant information of
 * the material quantity.
 * @param buffer_type interface description buffer type to be used (default:
 * Base).
 */
InterfaceFieldQuantity::InterfaceFieldQuantity(
    UnitHandler const &unit_handler, MaterialManager const &material_manager,
    std::string const &quantity_name, std::array<bool, 3> const output_flags,
    InterfaceFieldQuantityData const &quantity_data,
    InterfaceDescriptionBufferType const buffer_type)
    : // Start initializer list
      OutputQuantity(unit_handler, material_manager,
                     quantity_name +
                         SuffixOfInterfaceDescriptionBufferType(buffer_type),
                     output_flags, quantity_data.GetDimensions()),
      quantity_data_(quantity_data), buffer_type_(buffer_type) {
  /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief See base class definition.
 */
void InterfaceFieldQuantity::DoComputeCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter) const {

  // extract the correct field type for shorter access
  InterfaceFieldType const field_type = quantity_data_.field_type_;

  // Change calls for interface and non-interface nodes
  if (node.HasLevelset()) {
    // local counter
    unsigned long long int local_counter = 0;

    // Loop through all components
    for (std::size_t component = 0;
         component < quantity_data_.field_indices_.size(); component++) {
      // Get the correct field index
      unsigned int const field_index = quantity_data_.field_indices_[component];
      // Get buffers of the quantity
      double const(&field_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetFieldBuffer(field_type, field_index,
                                                  buffer_type_);
      // Dimensionalization factor for re-dimensionalization of variables
      double const dimensionalization_factor =
          unit_handler_.DimensionalizeValue(
              1.0, IF::FieldUnit(field_type, field_index));

      // set the local counter on original cell_data_counter + the component
      local_counter = cell_data_counter + component;

      // Loop through all internal cells of block
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
          for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
            cell_data[local_counter] =
                field_buffer[i][j][k] * dimensionalization_factor;

            local_counter += quantity_data_.field_indices_.size();
          }
        }
      }
    }
  } else {
    // If interface tags are included the signum is taken from the corner
    // interface tag
    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
    double const pre_factor = quantity_data_.use_interface_tags_for_default_
                                  ? Signum(interface_tags[0][0][0])
                                  : 1.0;
    double const default_value = pre_factor * quantity_data_.default_value_;

    // Add default value for non-interface blocks
    // Loop through all components
    for (std::size_t component = 0;
         component < quantity_data_.field_indices_.size(); component++) {
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
          for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
            cell_data[cell_data_counter++] = default_value;
          }
        }
      }
    }
  }
}

/**
 * @brief See base class definition.
 */
void InterfaceFieldQuantity::DoComputeDebugCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter, MaterialName const) const {

  // extract the correct field type, field index and unit for shorter access
  InterfaceFieldType const field_type = quantity_data_.field_type_;

  // Differ between interface and non-interface blocks
  if (node.HasLevelset()) {
    // local counter
    unsigned long long int local_counter = 0;

    // Loop through all components
    for (std::size_t component = 0;
         component < quantity_data_.field_indices_.size(); component++) {
      // Get the correct field index
      unsigned int const field_index = quantity_data_.field_indices_[component];
      // Get buffers of the quantity
      double const(&field_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetFieldBuffer(field_type, field_index,
                                                  buffer_type_);
      // Dimensionalization factor for re-dimensionalization of variables
      double const dimensionalization_factor =
          unit_handler_.DimensionalizeValue(
              1.0, IF::FieldUnit(field_type, field_index));

      // set the local counter on original cell_data_counter + the component
      local_counter = cell_data_counter + component;

      // Loop through all cells of block
      for (unsigned int k = 0; k < CC::TCZ(); ++k) {
        for (unsigned int j = 0; j < CC::TCY(); ++j) {
          for (unsigned int i = 0; i < CC::TCX(); ++i) {
            cell_data[local_counter] =
                field_buffer[i][j][k] * dimensionalization_factor;

            local_counter += quantity_data_.field_indices_.size();
          }
        }
      }
    }
  } else {
    // If interface tags are included the signum is taken from the corner
    // interface tag
    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

    // Add default value for non-interface blocks
    // Loop through all components
    for (std::size_t component = 0;
         component < quantity_data_.field_indices_.size(); component++) {
      for (unsigned int k = 0; k < CC::TCZ(); ++k) {
        for (unsigned int j = 0; j < CC::TCY(); ++j) {
          for (unsigned int i = 0; i < CC::TCX(); ++i) {
            double const pre_factor =
                quantity_data_.use_interface_tags_for_default_
                    ? double(interface_tags[i][j][k])
                    : 1.0;
            cell_data[cell_data_counter++] =
                pre_factor * quantity_data_.default_value_;
          }
        }
      }
    }
  }
}
