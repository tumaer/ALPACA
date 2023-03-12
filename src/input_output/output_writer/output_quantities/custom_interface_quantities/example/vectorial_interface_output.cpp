//===----------------- vectorial_interface_output.cpp ---------------------===//
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
#include "input_output/output_writer/output_quantities/custom_interface_quantities/example/vectorial_interface_output.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief constructor to create a vectorial interface output quantity.
 * @param unit_handler Unit handler class for dimensionalization.
 * @param material_manager Material manager for accessing material data.
 * @param quantity_name Name of the quantity which is displayed in the ParaView
 * cell-data list.
 * @param output_flags Flags for which output type an output is written (0:
 * standard, 1: interface, 2:debug).
 *
 * @note {row, colmun} = {DTI( CC::DIM() ),1} marks that the quantity is a
 * vector with components equal to current dimension of simulation, everything
 * else is treated as a matrix.
 * @note In ParaView a vector will be displayed by its name and four scalar
 * entries (X-, Y-, Z-Component and the total magnitude).
 */
VectorialInterfaceOutput::VectorialInterfaceOutput(
    UnitHandler const &unit_handler, MaterialManager const &material_manager,
    std::string const &quantity_name, std::array<bool, 3> const output_flags)
    : OutputQuantity(unit_handler, material_manager, quantity_name,
                     output_flags, {DTI(CC::DIM()), 1}) {
  /** Empty besides initializer list */
}

/**
 * @brief see base class definition.
 */
void VectorialInterfaceOutput::DoComputeCellData(
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
    // double (pre_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()][dimensions_[0]];

    // Loop through all internal cells to fill the vector appropriately
    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          for (unsigned int row = 0; row < dimensions_[0]; row++) {
            /**
             * Option 1: Assign the correct value by indexing
             */
            // cell_data[cell_data_counter++] = pre_buffer[i][j][k][row] *
            // dimensionalization_factor;

            /**
             * Option 2: No preparation has been done. Carry out computation
             * here
             */
            // Use data from the negative material buffer
            cell_data[cell_data_counter++] =
                2.0 * double(row + 1) * dimensionalization_factor;
          }

          /**
           * Option 3: Compute the vector for a single cell and assign it
           * properly
           */
          // std::array<double, dimensions[0]> vector = ComputeVector();
          // // Assign here the matrix tensor appropriately
          // cell_data[cell_data_counter] = vector[0];
          // cell_data[cell_data_counter + 1] = vector[1];
          // // ...
          // cell_data[cell_data_counter + dimensions[0] - 1] =
          // vector[dimensions[0] - 1];

          // // Increment the counter
          // cell_data_counter += dimensions[0];
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
          for (unsigned int row = 0; row < dimensions_[0]; row++) {
            cell_data[cell_data_counter++] =
                2.0 * double(row + 1) * dimensionalization_factor;
          }
        }
      }
    }
  }
}

/**
 * @brief see base class definition.
 */
void VectorialInterfaceOutput::DoComputeDebugCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter, MaterialName const) const {

  /**
   * Use the unit handler to specify the correct dimensionalization factor for
   * the quantity
   */
  double const dimensionalization_factor = 1.0;

  // Differ between nodes that contain the material to write the correct data,
  // otherwise write a default value
  if (node.HasLevelset()) {

    /**
     * Declare here all buffers that are used for this quantity from the
     * material fields. If gradients are used, remember to limit the
     * computations later by the interface tags to ensure that a gradient exists
     * or gives reasonable values.
     *
     * NOTE: Never ever change the loop structure, since the total number of
     * elements are required.
     *
     * See for options above in non-Debug mode
     */
    for (unsigned int k = 0; k < CC::TCZ(); ++k) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int i = 0; i < CC::TCX(); ++i) {
          for (unsigned int row = 0; row < dimensions_[0]; row++) {
            /**
             * Simply add to the final data vector
             */
            cell_data[cell_data_counter++] =
                2.0 * double(row + 1) * dimensionalization_factor;
          }
        }
      }
    }
  } else {
    // otherwise use debug_default_value
    for (unsigned int k = 0; k < CC::TCZ(); ++k) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int i = 0; i < CC::TCX(); ++i) {
          for (unsigned int row = 0; row < dimensions_[0]; row++) {
            /**
             * Default value if material is not contained in node
             */
            cell_data[cell_data_counter++] = -1.0;
          }
        }
      }
    }
  }
}
