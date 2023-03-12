//===------------------ instantiation_unit_handler.cpp --------------------===//
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
#include "instantiation/instantiation_unit_handler.h"
#include "input_output/input_reader.h"

namespace Instantiation {

/**
 * @brief Instantiates the complete unit handler class with the given input
 * reader.
 * @param input_reader Reader that provides access to the full data of the input
 * file.
 * @return The fully instantiated UnitHandler class.
 */
UnitHandler InstantiateUnitHandler(InputReader const &input_reader) {

  // read data
  double const reference_density =
      input_reader.GetDimensionalizationReader().ReadReferenceDensity();
  double const reference_velocity =
      input_reader.GetDimensionalizationReader().ReadReferenceVelocity();
  double const reference_length =
      input_reader.GetDimensionalizationReader().ReadReferenceLength();
  double const reference_temperature =
      input_reader.GetDimensionalizationReader().ReadReferenceTemperature();

  // Log data
  LogWriter &logger = LogWriter::Instance();
  logger.LogMessage(" ");
  logger.LogMessage("Dimensionalization parameter:");
  logger.LogMessage(
      StringOperations::Indent(2) + "Density reference    : " +
      StringOperations::ToScientificNotationString(reference_density, 9));
  logger.LogMessage(
      StringOperations::Indent(2) + "Velocity reference   : " +
      StringOperations::ToScientificNotationString(reference_velocity, 9));
  logger.LogMessage(
      StringOperations::Indent(2) + "Length reference     : " +
      StringOperations::ToScientificNotationString(reference_length, 9));
  logger.LogMessage(
      StringOperations::Indent(2) + "Temperature reference: " +
      StringOperations::ToScientificNotationString(reference_temperature, 9));
  logger.LogMessage(" ");

  // Initialize the unit handler
  return UnitHandler(reference_density, reference_velocity, reference_length,
                     reference_temperature);
}
} // namespace Instantiation
