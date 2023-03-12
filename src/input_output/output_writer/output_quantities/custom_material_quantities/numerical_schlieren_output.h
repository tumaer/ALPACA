//===----------------- numerical_schlieren_output.h -----------------------===//
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
#ifndef NUMERICAL_SCHLIEREN_OUTPUT_H
#define NUMERICAL_SCHLIEREN_OUTPUT_H

#include "input_output/output_writer/output_quantity.h"
#include "topology/node.h"

/**
 * @brief The NumericalSchlierenOutput class handles the output of the density
 * gradient into a output file (currently Xdmf + HDF5). NumericalSchlierenOutput
 * must not change any data.
 */
class NumericalSchlierenOutput : public OutputQuantity {

private:
  // Compute functions required from base class
  void
  DoComputeCellData(Node const &node, std::vector<double> &cell_data,
                    unsigned long long int &cell_data_counter) const override;
  void DoComputeDebugCellData(Node const &node, std::vector<double> &cell_data,
                              unsigned long long int &cell_data_counter,
                              MaterialName const material) const override;

public:
  NumericalSchlierenOutput() = delete;
  explicit NumericalSchlierenOutput(UnitHandler const &unit_handler,
                                    MaterialManager const &material_manager,
                                    std::string const &quantity_name,
                                    std::array<bool, 3> const output_flags);
  virtual ~NumericalSchlierenOutput() = default;
  NumericalSchlierenOutput(NumericalSchlierenOutput const &) = delete;
  NumericalSchlierenOutput &
  operator=(NumericalSchlierenOutput const &) = delete;
  NumericalSchlierenOutput(NumericalSchlierenOutput &&) = delete;
  NumericalSchlierenOutput &operator=(NumericalSchlierenOutput &&) = delete;
};

#endif // NUMERICAL_SCHLIEREN_OUTPUT_H
