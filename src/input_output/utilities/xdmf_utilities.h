//===------------------------- xdmf_utilities.h ---------------------------===//
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
#ifndef Xdmf_UTILITIES_H
#define Xdmf_UTILITIES_H

#include <array>
#include <string>

/**
 * @brief The XdmfUtilities serves as a helper class to provide the complete set
 * of properly formatted strings required for writing a xdmf file.
 */
namespace XdmfUtilities {

std::string TimeDataItem(double const output_time);
std::string DataItemString(std::string const &hdf5_filename,
                           std::string const &item_name,
                           unsigned long long int const number_of_cells,
                           std::array<unsigned int, 2> const &dimensions);
std::string TopologyString(std::string const &data_item,
                           unsigned long long int const number_of_cells);
std::string GeometryString(std::string const &data_item,
                           unsigned long long int const number_of_vertices);
std::string ScalarAttributeString(std::string const &attribute_name,
                                  std::string const &data_item);
std::string VectorAttributeString(std::string const &attribute_name,
                                  std::string const &data_item);
std::string MatrixAttributeString(std::string const &attribute_name,
                                  std::string const &data_item);
std::string TensorAttributeString(std::string const &attribute_name,
                                  std::string const &data_item);
std::string SpatialDataInformation(std::string const &spatial_data_name,
                                   std::string const &spatial_data_information);
std::string HeaderInformation(std::string const &data_name);
std::string FooterInformation();

} // namespace XdmfUtilities

#endif // Xdmf_UTILITIES_H
