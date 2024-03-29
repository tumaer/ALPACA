//===-------------- matrix_tensor_interface_output.cpp --------------------===//
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
#include "input_output/output_writer/output_quantities/custom_interface_quantities/example/matrix_tensor_interface_output.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief constructor to create a matrix/tensor interface output quantity.
 * @param unit_handler Unit handler class for dimensionalization.
 * @param material_manager Material manager for accessing material data.
 * @param quantity_name Name of the quantity which is displayed in the ParaView
 * cell-data list.
 * @param output_flags Flags for which output type an output is written (0:
 * standard, 1: interface, 2:debug).
 *
 * @note {row, colmun} = {DTI( CC::DIM(),1} marks that the quantity is a matrix
 * of current dimension and 2. For tensor { DTI( CC::DIM() ), DTI( CC::DIM() ) }
 * should be used.
 * @note In ParaView a matrix/tensor will be displayed by its name and m * n + 1
 * scalar entries (for each component + magnitude). The components are numbered
 * from 0 to m*n-1, where the ordering columns followed by rows (for a 2 x 2 =>
 * 0 : 00, 1: 01, 2: 10, 3: 11).
 */
MatrixTensorInterfaceOutput::MatrixTensorInterfaceOutput(
    UnitHandler const &unit_handler, MaterialManager const &material_manager,
    std::string const &quantity_name, std::array<bool, 3> const output_flags)
    : OutputQuantity(unit_handler, material_manager, quantity_name,
                     output_flags, {DTI(CC::DIM()), 2}) {
  /** Empty besides initializer list */
}

/**
 * @brief see base class definition.
 */
void MatrixTensorInterfaceOutput::DoComputeCellData(
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
    // double
    // (pre_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()][dimensions_[0]][dimensions_[1]];

    // Loop through all internal cells to fill the vector appropriately
    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          for (unsigned int row = 0; row < dimensions_[0]; row++) {
            for (unsigned int col = 0; col < dimensions_[1]; col++) {

              /**
               * Option 1: Assign the correct value by indexing
               */
              // cell_data[cell_data_counter++] = pre_buffer[i][j][k][row][col]
              // * dimensionalization_factor;

              /**
               * Option 2: No preparation has been done. Carry out computation
               * here
               */
              // Use data from the negative material buffer
              cell_data[cell_data_counter++] =
                  (2.0 * double(row + 1) + double(col + 1)) *
                  dimensionalization_factor;
            }
          }

          /**
           * Option 3: Compute the matrix/tensor for a single cell and assign it
           * properly
           */
          // std::array< std::array< double, dimensions[1]> dimensions[0]>
          // matrix_tensor = ComputeMatrixTensor();
          // // Assign here the matrix tensor appropriately
          // cell_data[cell_data_counter] = matrix_tensor[0][0];
          // cell_data[cell_data_counter + 1] = matrix_tensor[1][0];
          // // ...
          // cell_data[cell_data_counter + dimensions[0] * dimensions[1] - 1] =
          // matrix_tensor[dimensions[0] - 1][dimensions[1] - 1];

          // // Increment the counter
          // cell_data_counter += dimensions[0] * dimensions[1];
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
            for (unsigned int col = 0; col < dimensions_[1]; col++) {
              cell_data[cell_data_counter++] =
                  2.0 * (2.0 * double(row + 1) + double(col + 1)) *
                  dimensionalization_factor;
            }
          }
        }
      }
    }
  }
}

/**
 * @brief see base class definition.
 */
void MatrixTensorInterfaceOutput::DoComputeDebugCellData(
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
            for (unsigned int col = 0; col < dimensions_[1]; col++) {
              /**
               * Simply add to the final data vector
               */
              cell_data[cell_data_counter++] =
                  (2.0 * double(row + 1) + double(col + 1)) *
                  dimensionalization_factor;
            }
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
            for (unsigned int col = 0; col < dimensions_[1]; col++) {
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
}
