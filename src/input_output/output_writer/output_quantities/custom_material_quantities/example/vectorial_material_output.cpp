//===---------------- vectorial_material_output.cpp -----------------------===//
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
#include "input_output/output_writer/output_quantities/custom_material_quantities/example/vectorial_material_output.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief constructor to create a vectorial material output quantity.
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
VectorialMaterialOutput::VectorialMaterialOutput(
    UnitHandler const &unit_handler, MaterialManager const &material_manager,
    std::string const &quantity_name, std::array<bool, 3> const output_flags)
    : OutputQuantity(unit_handler, material_manager, quantity_name,
                     output_flags, {DTI(CC::DIM()), 1}) {
  /** Empty besides initializer list */
}

/**
 * @brief see base class definition.
 */
void VectorialMaterialOutput::DoComputeCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter) const {

  /**
   * Use the unit handler to specify the correct dimensionalization factor for
   * the quantity
   */
  double const dimensionalization_factor = 1.0;

  // Different behavior dependent on interface presence or not
  if (node.HasLevelset()) {
    // Use the interface tags and levelset to differ between different states
    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
    double const(&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetBaseBuffer(InterfaceDescription::Levelset);

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
          // Loop through all rows and columns
          for (unsigned int row = 0; row < dimensions_[0]; row++) {

            /**
             * Option 1: Assign the correct value by indexing
             */
            // cell_data[cell_data_counter++] = pre_buffer[i][j][k][row] *
            // dimensionalization_factor;

            /**
             * Option 2: No preparation has been done. Differ between materials
             * at this stage
             */

            // Use data from the negative material buffer
            if (interface_tags[i][j][k] < 0 ||
                (std::abs(interface_tags[i][j][k]) <= ITTI(IT::NewCutCell) &&
                 levelset[i][j][k] < 0.0)) {
              cell_data[cell_data_counter++] =
                  2.0 * double(row + 1) * dimensionalization_factor;
            }
            // otherwise positive
            else {
              cell_data[cell_data_counter++] =
                  -2.0 * double(row + 1) * dimensionalization_factor;
            }
          }

          /**
           * Option 3: Compute the matrix/tensor for a single cell and assign it
           * properly
           */
          // std::array< std::array< double, dimensions[1]> dimensions[0]>
          // vector = ComputeVector();
          // // Assign here the matrix tensor appropriately
          // cell_data[cell_data_counter] = vector[0];
          // cell_data[cell_data_counter + 1] = vector[1];
          // // ...
          // cell_data[cell_data_counter + dimensions[0] - 1] =
          // vector[dimensions[0] - 1];

          // // Increment the counter
          // cell_data_counter += dimensions[0] * dimensions[1];
        }
      }
    }
  } else {
    // Extract the single material present on this node
    // No interface node -> interface tags/material is the same everywhere
    MaterialName const material = node.GetSinglePhaseMaterial();
    Block const &block = node.GetPhaseByMaterial(material);
    (void)block;
    /**
     * Declare here all buffers that are used for this quantity from the
     * material fields
     */
    // Fill the data vector with the appropriate values
    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          /**
           * Add data to vector (see above)
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
 *
 * @note Attention: In case prime state, parameter  variables are used, pay
 * attention that they only exist on leave nodes. In case a division is made on
 * non-leave nodes a floating point exception is caused. Therefore, only use the
 * debug output if it is ensured that this cannot happen. Conservatives can be
 * used since they are present on all nodes.
 */
void VectorialMaterialOutput::DoComputeDebugCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter,
    MaterialName const material) const {

  /**
   * Use the unit handler to specify the correct dimensionalization factor for
   * the quantity
   */
  double const dimensionalization_factor = 1.0;

  // Differ between nodes that contain the material to write the correct data,
  // otherwise write a default value
  if (node.ContainsMaterial(material)) {

    /**
     * Declare here all buffers that are used for this quantity from the
     * material fields. If gradients are used, remember to limit the
     * computations later by the interface tags to ensure that a gradient exists
     * or gives reasonable values.
     *
     * Never ever change the loop structure, since the total number of elements
     * are required.
     *
     * See details above in non-Debug mode
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
