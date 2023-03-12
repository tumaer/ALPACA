//===------------------ scalar_interface_output.cpp -----------------------===//
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
#include "input_output/output_writer/output_quantities/custom_interface_quantities/example/scalar_interface_output.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief constructor to create a scalar interface output quantity.
 * @param unit_handler Unit handler class for dimensionalization.
 * @param material_manager Material manager for accessing material data.
 * @param quantity_name Name of the quantity which is displayed in the ParaView
 * cell-data list.
 * @param output_flags Flags for which output type an output is written (0:
 * standard, 1: interface, 2:debug).
 *
 * @note {row, colmun} = {1,1} marks that the quantity is a scalar.
 */
ScalarInterfaceOutput::ScalarInterfaceOutput(
    UnitHandler const &unit_handler, MaterialManager const &material_manager,
    std::string const &quantity_name, std::array<bool, 3> const output_flags)
    : OutputQuantity(unit_handler, material_manager, quantity_name,
                     output_flags, {1, 1}) {
  /** Empty besides initializer list */
}

/**
 * @brief see base class definition.
 */
void ScalarInterfaceOutput::DoComputeCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter) const {

  /**
   * Use the unit handler to specify the correct dimensionalization factor for
   * the quantity
   */
  double const dimensionalization_factor = 1.0;

  // Different behavior dependent on interface presence or not
  if (node.HasLevelset()) {

    /**
     * Declare here all buffers that are used for this quantity from the
     * material fields
     */
    // double const (&positive_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()];
    // double const (&negative_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()];

    /**
     * Option 1: Here, carry out pre-operations to fill a complete new buffer
     * (e.g. one that contains only the real-material properties depending on
     * the interface tags, which could be required for the computation of
     * derivatives).
     */
    // double (pre_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()];

    // Loop through all internal cells to fill the vector appropriately
    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          /**
           * Option 1: Assign the correct value by indexing
           */
          // cell_data[cell_data_counter++] = pre_buffer[i][j][k] *
          // dimensionalization_factor;

          /**
           * Option 2: No preparation has been done. Carry out computation here
           */
          // Use data from the negative material buffer
          cell_data[cell_data_counter++] = 1.0 * dimensionalization_factor;
        }
      }
    }
  } else {
    /**
     * Declare here all buffers that are used for this quantity from the
     * material fields
     */

    // Fill the data vector with the appropriate values
    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          /**
           * Simply add to the final data vector
           */
          cell_data[cell_data_counter++] = -1.0 * dimensionalization_factor;
        }
      }
    }
  }
}

/**
 * @brief see base class definition.
 */
void ScalarInterfaceOutput::DoComputeDebugCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter, MaterialName const) const {

  /**
   * Use the unit handler to specify the correct dimensionalization factor for
   * the quantity
   */
  double const dimensionalization_factor = 1.0;

  // Change calls for interface and non-interface nodes
  if (node.HasLevelset()) {

    /**
     * Declare here all buffers that are used for this quantity from the
     * interface fields. If gradients are used Remember to limit the output to
     * internal cells.
     *
     * Never ever change the loop structure, since the total number of elements
     * are required.
     *
     * See for options above in non-Debug mode
     */

    for (unsigned int k = 0; k <= CC::TCZ(); ++k) {
      for (unsigned int j = 0; j <= CC::TCY(); ++j) {
        for (unsigned int i = 0; i <= CC::TCX(); ++i) {
          /**
           * Add the value to the data vector (with further operations)
           */
          cell_data[cell_data_counter++] = 1.0 * dimensionalization_factor;
        }
      }
    }

  } else {
    // Add default value for non-interface blocks
    for (unsigned int k = 0; k <= CC::TCZ(); ++k) {
      for (unsigned int j = 0; j <= CC::TCY(); ++j) {
        for (unsigned int i = 0; i <= CC::TCX(); ++i) {
          /**
           * Add the default value to the output data
           */
          cell_data[cell_data_counter++] = -1.0;
        }
      }
    }
  }
}
