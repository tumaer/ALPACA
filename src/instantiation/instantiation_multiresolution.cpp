//===-------------- instantiation_multiresolution.cpp ---------------------===//
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
#include "instantiation/instantiation_multiresolution.h"

namespace Instantiation {

/**
 * @brief Creates a Multiresolution object.
 * @param input_reader Reader that provides access to the full data of the input
 * file.
 * @param topology_manager Class providing global (on all ranks) node
 * information.
 * @return A multiresolution object with thresholding conditions accroding to
 * given input.
 */
Multiresolution
InstantiateMultiresolution(InputReader const &input_reader,
                           TopologyManager const &topology_manager) {

  // Get the data local (for logging)
  unsigned int const maximum_level = topology_manager.GetMaximumLevel();
  unsigned int const epsilon_reference_level =
      input_reader.GetMultiResolutionReader().ReadEpsilonLevelReference();
  double const epsilon_reference =
      input_reader.GetMultiResolutionReader().ReadEpsilonReference();

  // Create the thresholder
  Thresholder thresholder =
      Thresholder(maximum_level, epsilon_reference_level, epsilon_reference);

  // Log data
  LogWriter &logger = LogWriter::Instance();
  logger.LogMessage(" ");
  if (maximum_level > 1) {
    // tmp string for maximum size determination
    std::string const level_string(std::to_string(maximum_level - 1) +
                                   std::to_string(maximum_level));

    logger.LogMessage("Epsilon Level 0 and Level 1" +
                      std::string(level_string.size() - 2, ' ') + " : " +
                      StringOperations::ToScientificNotationString(
                          thresholder.ThresholdOnLevel(1), 9));
    logger.LogMessage("Epsilon Level " + std::to_string(maximum_level - 1) +
                      " and Level " + std::to_string(maximum_level) + " : " +
                      StringOperations::ToScientificNotationString(
                          thresholder.ThresholdOnLevel(maximum_level), 9));
  } else if (maximum_level == 1) {
    logger.LogMessage("Homogenous Mesh - Level 1 cannot be coarsened. Provided "
                      "Epsilon is ignored");
  } else {
    logger.LogMessage("Homogenous Mesh - Provided Epsilon is ignored ");
  }
  logger.LogMessage(" ");

  // Return the created multiresolution with the thresholder
  return Multiresolution(std::move(thresholder));
}
} // namespace Instantiation
