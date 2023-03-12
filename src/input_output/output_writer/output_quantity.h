//===------------------------ output_quantity.h ---------------------------===//
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
#ifndef OUTPUT_QUANTITY_H
#define OUTPUT_QUANTITY_H

#include "input_output/output_writer/output_definitions.h"
#include "materials/material_definitions.h"
#include "materials/material_manager.h"
#include "topology/node.h"
#include "unit_handler.h"
#include <hdf5.h>

/**
 * @brief The OutputQuantity class handles the output to the filesystem in
 * Xdmf+HDF5 file format for ParaView. OutputQuantity must not change any data.
 * @note This class works as a base class where all field output quantities must
 * inherit from that write cell data to the appropriate file.
 */
class OutputQuantity {

protected:
  // unit_handler for the dimensionalization of variables
  UnitHandler const &unit_handler_;
  // material manager to allow computations for different materials and use
  // properties/eos data or pairing data
  MaterialManager const &material_manager_;
  // flags to specify where the quantities should be used (0: standard, 1:
  // interface, 2:debug)
  std::array<bool, 3> const output_flags_;
  // name of the complete quantity used in the hdf5 + xdmf file for naming
  std::string const quantity_name_;
  // dimensions of the quantity ( { 1 } : Scalar, { n } : vector, { n x m } :
  // matrix )
  std::array<unsigned int, 2> const dimensions_;

  /**
   * @brief Compute values for the data vector that is written to the hdf5 file
   * (standard, interface mode).
   * @param node Node data that should be added.
   * @param cell_data_counter Actual position of the index in the cell data
   * vector.
   * @param cell_data Vector in which the data is written.
   *
   * @note pure virtual function that MUST be implemented by all derived
   * classes.
   */
  virtual void
  DoComputeCellData(Node const &node, std::vector<double> &cell_data,
                    unsigned long long int &cell_data_counter) const = 0;

  /**
   * @brief Compute values for the data vector that is written to the hdf5 file
   * (debug mode).
   * @param node Node data that should be added.
   * @param cell_data_counter Actual position of the index in the cell data
   * vector.
   * @param cell_data Vector in which the data is written.
   * @param material Material that should be considered (only used for material
   * output quantities).
   *
   * @note pure virtual function that MUST be implemented by all derived classes
   */
  virtual void DoComputeDebugCellData(Node const &node,
                                      std::vector<double> &cell_data,
                                      unsigned long long int &cell_data_counter,
                                      MaterialName const material) const = 0;

  // protected constructor to prevent direct calls from outside
  explicit OutputQuantity(UnitHandler const &unit_handler,
                          MaterialManager const &material_manager,
                          std::string const &quantity_name,
                          std::array<bool, 3> const &output_flags,
                          std::array<unsigned int, 2> const &dimensions);

public:
  OutputQuantity() = delete;
  virtual ~OutputQuantity() = default;
  OutputQuantity(OutputQuantity const &) = delete;
  OutputQuantity &operator=(OutputQuantity const &) = delete;
  OutputQuantity(OutputQuantity &&) = delete;
  OutputQuantity &operator=(OutputQuantity &&) = delete;

  // Functions that fill the complete cell data vector with appropriate values
  void
  ComputeCellData(std::vector<std::reference_wrapper<Node const>> const &nodes,
                  std::vector<double> &cell_data) const;
  void ComputeDebugCellData(
      std::vector<std::reference_wrapper<Node const>> const &nodes,
      std::vector<double> &cell_data,
      MaterialName const material = MaterialName::MaterialOne) const;

  // Creates and attribute string compliant with a xdmf reader
  std::string GetXdmfAttributeString(std::string const &filename,
                                     std::string const &group_name,
                                     hsize_t const number_of_values,
                                     std::string const prefix = "") const;

  // Additional return functions to provide data to outside
  bool IsActive(OutputType const output_type) const;
  std::string GetName() const;
  std::array<unsigned int, 2> GetDimensions() const;
};

#endif // OUTPUT_QUANTITY_H
