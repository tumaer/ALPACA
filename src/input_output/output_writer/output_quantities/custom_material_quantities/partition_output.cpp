//===---------------------- partition_output.cpp --------------------------===//
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
#include "input_output/output_writer/output_quantities/custom_material_quantities/partition_output.h"
#include "communication/mpi_utilities.h"

/**
 * @brief constructor to output the rank number.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @param material_manager Instance to access all material data.
 * @param quantity_name Name of the quantity that is displayed in the ParaView
 * cell data list.
 * @param output_flags Flags of the output type that is written (0: standard, 1:
 * interface, 2:debug).
 *
 * @note {row, colmun} = {1,1} marks that the quantity is a scalar.
 */
PartitionOutput::PartitionOutput(UnitHandler const &unit_handler,
                                 MaterialManager const &material_manager,
                                 std::string const &quantity_name,
                                 std::array<bool, 3> const output_flags)
    : OutputQuantity(unit_handler, material_manager, quantity_name,
                     output_flags, {1, 1}),
      rank_in_double_format_(double(MpiUtilities::MyRankId())) {
  /** Empty besides initializer list */
}

/**
 * @brief see base class definition.
 */
void PartitionOutput::DoComputeCellData(
    Node const &, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter) const {

  // Loop through number of internal cells
  for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
        cell_data[cell_data_counter++] = rank_in_double_format_;
      }
    }
  }
}

/**
 * @brief see base class definition.
 */
void PartitionOutput::DoComputeDebugCellData(
    Node const &, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter, MaterialName const) const {

  /** Assign the correct rank to the data vector */
  for (unsigned int k = 0; k < CC::TCZ(); ++k) {
    for (unsigned int j = 0; j < CC::TCY(); ++j) {
      for (unsigned int i = 0; i < CC::TCX(); ++i) {
        cell_data[cell_data_counter++] = rank_in_double_format_;
      }
    }
  }
}
