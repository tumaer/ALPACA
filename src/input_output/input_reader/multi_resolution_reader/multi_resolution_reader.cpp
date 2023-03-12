//===----------------- multi_resolution_reader.cpp ------------------------===//
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
#include "input_output/input_reader/multi_resolution_reader/multi_resolution_reader.h"

#include <numeric>

#include "enums/dimension_definition.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Gives the checked number of nodes on level zero for a given direction.
 * @param direction Direction identifier for which the number of nodes should be
 * read.
 * @return number of nodes.
 */
unsigned int
MultiResolutionReader::ReadNumberOfNodes(Direction const direction) const {
  // Read the block number and check on consistency
  int const number_of_nodes(DoReadNumberOfNodes(direction));
  if (number_of_nodes <= 0) {
    throw std::invalid_argument(
        "At least one block must be present on level zero in each direction!");
  }
  if (number_of_nodes > 128) {
    throw std::invalid_argument(
        "Block number on level zero must not exceed 128!");
  }
  return static_cast<unsigned int>(number_of_nodes);
}

/**
 * @brief Gives the checked size of a node on level zero.
 * @return node size on level zero.
 */
double MultiResolutionReader::ReadNodeSizeOnLevelZero() const {
  // Read the node size value and check on consistency
  double const node_size(DoReadNodeSizeOnLevelZero());
  if (node_size <= 0.0) {
    throw std::invalid_argument(
        "Node size on level zero must be greater zero!");
  }
  return node_size;
}

/**
 * @brief Gives the checked maximum level used for the simulation.
 * @return maximum level of the simulation.
 */
unsigned int MultiResolutionReader::ReadMaximumLevel() const {
  // Read the value and check on consistency
  int const maximum_level(DoReadMaximumLevel());
  if (maximum_level > static_cast<int>(CC::AMNL())) {
    throw std::invalid_argument("Maximum level must NOT be larger than " +
                                std::to_string(CC::AMNL()) + "!");
  }
  if (maximum_level < 0) {
    throw std::invalid_argument("Maximum level must NOT be below zero!");
  }
  return static_cast<unsigned int>(maximum_level);
}

/**
 * @brief Gives the checked reference epsilon value used for the refinement
 * criterion.
 * @return epsilon refernce value.
 */
double MultiResolutionReader::ReadEpsilonReference() const {
  // Read the reference value and check on consistency
  double const reference(DoReadEpsilonReference());
  if (reference <= 0.0) {
    throw std::invalid_argument("Epsilon reference must be larger than zero!");
  }
  return reference;
}

/**
 * @brief Gives the checked reference level where the epsilon reference is
 * enforced.
 * @return level of reference.
 */
unsigned int MultiResolutionReader::ReadEpsilonLevelReference() const {
  // Read the reference level and check on consistency
  int const reference_level(DoReadEpsilonLevelReference());
  if (reference_level < 0) {
    throw std::invalid_argument(
        "Level of epsilon reference must NOT be below zero!");
  }
  return static_cast<unsigned int>(reference_level);
}
