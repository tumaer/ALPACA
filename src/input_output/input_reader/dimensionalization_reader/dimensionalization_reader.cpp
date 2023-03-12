//===------------------ dimensionalization_reader.cpp ---------------------===//
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
#include "input_output/input_reader/dimensionalization_reader/dimensionalization_reader.h"
#include <stdexcept>

/**
 * @brief Gives the checked reference length of the dimensionalization
 * parameter.
 * @return length refernce value.
 */
double DimensionalizationReader::ReadReferenceLength() const {
  // Get the reference value and make consistency check
  double const reference_value(DoReadReferenceLength());
  if (reference_value <= 0.0) {
    throw std::invalid_argument("Reference length must be larger than zero!");
  }
  return reference_value;
}

/**
 * @brief Gives the checked reference velocity of the dimensionalization
 * parameter.
 * @return velocity refernce value.
 */
double DimensionalizationReader::ReadReferenceVelocity() const {
  // Get the reference value and make consistency check
  double const reference_value(DoReadReferenceVelocity());
  if (reference_value <= 0.0) {
    throw std::invalid_argument("Reference velocity must be larger than zero!");
  }
  return reference_value;
}

/**
 * @brief Gives the checked reference density of the dimensionalization
 * parameter.
 * @return density refernce value.
 */
double DimensionalizationReader::ReadReferenceDensity() const {
  // Get the reference value and make consistency check
  double const reference_value(DoReadReferenceDensity());
  if (reference_value <= 0.0) {
    throw std::invalid_argument("Reference density must be larger than zero!");
  }
  return reference_value;
}

/**
 * @brief Gives the checked reference temperature of the dimensionalization
 * parameter.
 * @return temperature refernce value.
 */
double DimensionalizationReader::ReadReferenceTemperature() const {
  // Get the reference value and make consistency check
  double const reference_value(DoReadReferenceTemperature());
  if (reference_value <= 0.0) {
    throw std::invalid_argument(
        "Reference temperature must be larger than zero!");
  }
  return reference_value;
}
