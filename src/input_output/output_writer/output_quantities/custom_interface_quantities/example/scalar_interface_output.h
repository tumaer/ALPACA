//===------------------- scalar_interface_output.h ------------------------===//
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
#ifndef SCALAR_INTERFACE_OUTPUT_H
#define SCALAR_INTERFACE_OUTPUT_H

#include "input_output/output_writer/output_quantity.h"
#include "topology/node.h"

/**
 * @brief Example class for a scalar interface quantity that writes interface
 * data to the output. Interface data is defined as all data operating on a
 * interface block. The difference to material outputs is that those operate on
 * Material blocks.
 *
 * To define a new scalar interface quantity this file can be used to define the
 * appropriate functions. Simply take the current loop structure and add
 * required operations at the already specified locations. To use the output
 * quantity the following two steps are needed:
 *    1. An output setting must be written in the file
 * "user_specifications/output_constants.h" (see examples there)
 *    2. Add the constructor call in the "instantiation_output_writer.cpp"
 * function.
 */
class ScalarInterfaceOutput : public OutputQuantity {

private:
  // Compute functions required from base class
  void
  DoComputeCellData(Node const &node, std::vector<double> &cell_data,
                    unsigned long long int &cell_data_counter) const override;
  void DoComputeDebugCellData(Node const &node, std::vector<double> &cell_data,
                              unsigned long long int &cell_data_counter,
                              MaterialName const material) const override;

public:
  ScalarInterfaceOutput() = delete;
  explicit ScalarInterfaceOutput(UnitHandler const &unit_handler,
                                 MaterialManager const &material_manager,
                                 std::string const &quantity_name,
                                 std::array<bool, 3> const output_flags);
  virtual ~ScalarInterfaceOutput() = default;
  ScalarInterfaceOutput(ScalarInterfaceOutput const &) = delete;
  ScalarInterfaceOutput &operator=(ScalarInterfaceOutput const &) = delete;
  ScalarInterfaceOutput(ScalarInterfaceOutput &&) = delete;
  ScalarInterfaceOutput &operator=(ScalarInterfaceOutput &&) = delete;
};

#endif // SCALAR_INTERFACE_OUTPUT_H
