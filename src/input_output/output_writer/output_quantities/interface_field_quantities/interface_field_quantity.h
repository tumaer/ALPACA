//===------------------ interface_field_quantity.h ------------------------===//
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
#ifndef INTERFACE_FIELD_QUANTITY_H
#define INTERFACE_FIELD_QUANTITY_H

#include "input_output/output_writer/output_quantities/interface_field_quantities/interface_field_quantity_definitions.h"
#include "input_output/output_writer/output_quantity.h"
#include "unit_handler.h"

/**
 * @brief The InterfaceFieldQuantity class handles the output to the filesystem
 * in Xdmf+HDF5 file format for ParaView. InterfaceFieldQuantity must not change
 * any data.
 *
 *        This quantity is used to write all interface fields to the hdf5/xdmf
 * file if desired. For activating different fields use the "output_constants.h"
 * file. If a new interface field needs to be added, refer to the
 * "interface_field_quantity_definitions.h". Here, nothing has to be changed.
 */
class InterfaceFieldQuantity : public OutputQuantity {

private:
  // struct containing all data required for the output
  InterfaceFieldQuantityData const quantity_data_;
  // interface description buffer type
  InterfaceDescriptionBufferType const buffer_type_;

  // Compute functions required from base class
  void
  DoComputeCellData(Node const &node, std::vector<double> &cell_data,
                    unsigned long long int &cell_data_counter) const override;
  void DoComputeDebugCellData(Node const &node, std::vector<double> &cell_data,
                              unsigned long long int &cell_data_counter,
                              MaterialName const material) const override;

public:
  InterfaceFieldQuantity() = delete;
  explicit InterfaceFieldQuantity(
      UnitHandler const &unit_handler, MaterialManager const &material_manager,
      std::string const &quantity_name, std::array<bool, 3> const output_flags,
      InterfaceFieldQuantityData const &quantity_data,
      InterfaceDescriptionBufferType const buffer_type =
          InterfaceDescriptionBufferType::Base);
  virtual ~InterfaceFieldQuantity() = default;
  InterfaceFieldQuantity(InterfaceFieldQuantity const &) = delete;
  InterfaceFieldQuantity &operator=(InterfaceFieldQuantity const &) = delete;
  InterfaceFieldQuantity(InterfaceFieldQuantity &&) = delete;
  InterfaceFieldQuantity &operator=(InterfaceFieldQuantity &&) = delete;
};

#endif // INTERFACE_FIELD_QUANTITY_H
