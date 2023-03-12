//===--------------------- instantiation_tree.cpp -------------------------===//
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
#include "instantiation/topology/instantiation_tree.h"

namespace Instantiation {

/**
 * @brief Instantiates the complete tree class with the given input classes.
 * @param input_reader Reader that provides access to the full data of the input
 * file.
 * @param topology_manager Class providing global (on all ranks) node
 * information.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @return The fully instantiated tree class.
 */
Tree InstantiateTree(InputReader const &input_reader,
                     TopologyManager &topology_manager,
                     UnitHandler const &unit_handler) {

  // Get the data for initializing and/or logging
  unsigned int const maximum_level = topology_manager.GetMaximumLevel();
  std::array<unsigned int, 3> const number_of_nodes =
      topology_manager.GetNumberOfNodesOnLevelZero();
  double const node_size_on_level_zero =
      input_reader.GetMultiResolutionReader().ReadNodeSizeOnLevelZero();

  // Logging of data
  LogWriter &logger = LogWriter::Instance();
  std::string tmp_string;

  logger.LogMessage(" ");
  // Domain size
  tmp_string = "Domain size                : " +
               StringOperations::ToScientificNotationString(
                   number_of_nodes[0] * node_size_on_level_zero, 6);
  tmp_string +=
      CC::DIM() != Dimension::One
          ? " x " + StringOperations::ToScientificNotationString(
                        number_of_nodes[1] * node_size_on_level_zero, 6)
          : "";
  tmp_string +=
      CC::DIM() == Dimension::Three
          ? " x " + StringOperations::ToScientificNotationString(
                        number_of_nodes[2] * node_size_on_level_zero, 6)
          : "";
  logger.LogMessage(tmp_string);
  // Internal cells per block
  logger.LogMessage("Internal cell per block    : " +
                    std::to_string(CC::ICX()));
  // Maximum level
  logger.LogMessage("Maximum level              : " +
                    std::to_string(maximum_level));
  // Resolution level zero
  tmp_string = "Resolution on level zero   : " +
               std::to_string(number_of_nodes[0] * CC::ICX());
  tmp_string += CC::DIM() != Dimension::One
                    ? " x " + std::to_string(number_of_nodes[1] * CC::ICX())
                    : "";
  tmp_string += CC::DIM() == Dimension::Three
                    ? " x " + std::to_string(number_of_nodes[2] * CC::ICX())
                    : "";
  logger.LogMessage(tmp_string + " internal cells");
  // Cell size on level zero
  logger.LogMessage("Cell size on level zero    : " +
                    StringOperations::ToScientificNotationString(
                        node_size_on_level_zero / double(CC::ICX()), 9));
  // Resolution level maximum
  tmp_string =
      "Resolution on maximum level: " +
      std::to_string(number_of_nodes[0] * CC::ICX() * (1 << maximum_level));
  tmp_string += CC::DIM() != Dimension::One
                    ? " x " + std::to_string(number_of_nodes[1] * CC::ICX() *
                                             (1 << maximum_level))
                    : "";
  tmp_string += CC::DIM() == Dimension::Three
                    ? " x " + std::to_string(number_of_nodes[2] * CC::ICX() *
                                             (1 << maximum_level))
                    : "";
  logger.LogMessage(tmp_string + " internal cells");
  // Cell size on level maximum
  logger.LogMessage("Cell size on maximum level : " +
                    StringOperations::ToScientificNotationString(
                        node_size_on_level_zero / double(CC::ICX()) /
                            double(1 << maximum_level),
                        9));
  logger.LogMessage(" ");

  // return the initialized tree
  return Tree(topology_manager, topology_manager.GetMaximumLevel(),
              unit_handler.NonDimensionalizeValue(node_size_on_level_zero,
                                                  UnitType::Length));
}
} // namespace Instantiation
