//===---------------------- output_quantity.cpp ---------------------------===//
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
#include "output_quantity.h"
#include "input_output/utilities/xdmf_utilities.h"

/**
 * @brief Explicit constructor to be used to create an output quantity.
 * @param unit_handler Handler to allow dimensionalization.
 * @param material_manager Manager for handling of the materials.
 * @param output_flags Flags specifying for which output type the quantity
 * should be written.
 * @param quantity_name Name of the given quantity used in the hdf5 and xdmf
 * file.
 * @param dimensions Dimensions of the given quantity ( {1,1} : scalar, {3,1}:
 * vector, {n,m}: matrix/tensor).
 */
OutputQuantity::OutputQuantity(UnitHandler const &unit_handler,
                               MaterialManager const &material_manager,
                               std::string const &quantity_name,
                               std::array<bool, 3> const &output_flags,
                               std::array<unsigned int, 2> const &dimensions)
    : // Start initializer list
      unit_handler_(unit_handler), material_manager_(material_manager),
      output_flags_(output_flags), quantity_name_(quantity_name),
      dimensions_(dimensions) {
  /** Empty constructor besides initializer list */
}

/**
 * @brief Gives the name of the quantity used in the xdmf and hdf5 file.
 * @return Name of the quantity.
 */
std::string OutputQuantity::GetName() const { return quantity_name_; }

/**
 * @brief Gives the dimensions of the quantity.
 * @return DImensions of the quantity.
 */
std::array<unsigned int, 2> OutputQuantity::GetDimensions() const {
  return dimensions_;
}

/**
 * @brief Checks whether the output quantity should be written for a given
 * output type (0: standard, 1: interface, 2: debug).
 * @param output_type Output type to be checked.
 * @return true if output should be written, false if not.
 */
bool OutputQuantity::IsActive(OutputType const output_type) const {
  return output_flags_[OTTI(output_type)];
}

/**
 * @brief Gives the appropriate attribute string for the xdmf file for this
 * quantity (differentiation between multidimensional, vectorial and scalar
 * quantities).
 * @param hdf5_filename HDF5 filename (without path) where the actual data has
 * been written to.
 * @param group_name Name of the group the data was written.
 * @param number_of_global_cells Number of cells used for the complete quantity
 * (globally on all ranks).
 * @param prefix Additional prefix that can be used for the attribute_name.
 * @return Compete Xdmf attribute string.
 */
std::string OutputQuantity::GetXdmfAttributeString(
    std::string const &hdf5_filename, std::string const &group_name,
    hsize_t const number_of_global_cells, std::string const prefix) const {
  // Get the data item
  std::string const data_item(XdmfUtilities::DataItemString(
      hdf5_filename, group_name + "/" + prefix + quantity_name_,
      number_of_global_cells, dimensions_));
  if (dimensions_.back() > 1 ||
      dimensions_.front() >
          DTI(CC::DIM())) { // multidimensional (second component dimension
                            // larger than one or first larger than the current
                            // dimension)
    // Both components are equal (nxn) and smaller or equal the dimension ->
    // Assumption that it is a tensor This is not relevant for the reading, but
    // allows sometimes special computations in ParaView
    if (dimensions_.back() == dimensions_.front() &&
        dimensions_.back() <= DTI(CC::DIM())) {
      return XdmfUtilities::TensorAttributeString(prefix + quantity_name_,
                                                  data_item);
    } else {
      return XdmfUtilities::MatrixAttributeString(prefix + quantity_name_,
                                                  data_item);
    }
  } else if (dimensions_.front() >
             1) { // vectorial (first component dimension larger than one and
                  // implicitly smaller than current dimension)
    return XdmfUtilities::VectorAttributeString(prefix + quantity_name_,
                                                data_item);
  } else { // scalar (both dimensions are one)
    return XdmfUtilities::ScalarAttributeString(prefix + quantity_name_,
                                                data_item);
  }
}

/**
 * @brief Computes the data of all internal cells and writes them into the
 * output data vector.
 * @param nodes Vector of nodes for which the output should be done.
 * @param cell_data Vector in which the full set of cell data is written.
 *
 * @note The cell_data vector stores its values in a pre-defined order. For a n
 * x m output quantity with and 3 cells the order is a follows: cell 1 | cell 2
 * |              cell 3                   | 00 01 ... 0m 10 ... 1m ... n0 ...
 * nm  | 00 01 ... 0m 10 ... 1m ... n0 ... nm  | 00 01 ... 0m 10 ... 1m ... n0
 * ... nm  |
 */
void OutputQuantity::ComputeCellData(
    std::vector<std::reference_wrapper<Node const>> const &nodes,
    std::vector<double> &cell_data) const {

  // Counter for the cell data vector
  unsigned long long int cell_data_counter = 0;
  // node counter
  unsigned int node_counter = 0;
  // factor that specifies the number of cell data entries that must be/are
  // written per node
  unsigned long long int cell_data_per_node =
      dimensions_[0] * dimensions_[1] * CC::ICX() * CC::ICY() * CC::ICZ();

  // Loop through all nodes to append correct number of data
  for (Node const &node : nodes) {
    // Ensures that the data counter is at the correct position for this node
    cell_data_counter = cell_data_per_node * node_counter;
    // Call the function that is implemented by the derived class
    DoComputeCellData(node, cell_data, cell_data_counter);
    // Increment the node counter
    node_counter++;
  }
}

/**
 * @brief Computes the data of all cells (internal + halos) and writes them into
 * the output data vector.
 * @param nodes Vector of nodes for which the output should be done.
 * @param cell_data Vector in which the full set of cell data is written.
 * @param material Material used for the output (not always used, only for
 * material output quantities).
 *
 * @note The cell_data vector stores its values in a pre-defined order. For a n
 * x m output quantity with and 3 cells the order is a follows: cell 1 | cell 2
 * |              cell 3                   | 00 01 ... 0m 10 ... 1m ... n0 ...
 * nm  | 00 01 ... 0m 10 ... 1m ... n0 ... nm  | 00 01 ... 0m 10 ... 1m ... n0
 * ... nm  |
 */
void OutputQuantity::ComputeDebugCellData(
    std::vector<std::reference_wrapper<Node const>> const &nodes,
    std::vector<double> &cell_data, MaterialName const material) const {

  // Counter for the cell data vector
  unsigned long long int cell_data_counter = 0;
  // node counter
  unsigned int node_counter = 0;
  // factor that specifies the number of cell data entries that must be/are
  // written per node
  unsigned long long int cell_data_per_node =
      dimensions_[0] * dimensions_[1] * CC::TCX() * CC::TCY() * CC::TCZ();

  // Loop through all nodes to append correct number of data
  for (Node const &node : nodes) {
    // Ensures that the data counter is at the correct position for this node
    cell_data_counter = cell_data_per_node * node_counter;
    // Call the function that is implemented by the derived class
    DoComputeDebugCellData(node, cell_data, cell_data_counter, material);
    // Increment the node counter
    node_counter++;
  }
}
