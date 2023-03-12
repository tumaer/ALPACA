//===------------------------ output_reader.cpp ---------------------------===//
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
#include "input_output/input_reader/output_reader/output_reader.h"

#include <algorithm>

/**
 * @brief Gives the checked time naming factor used for naming the output files.
 * @return factor used for the naming.
 */
double OutputReader::ReadTimeNamingFactor() const {
  // Read the time interval and check on consistency
  double const factor(DoReadTimeNamingFactor());
  if (factor <= 0.0) {
    throw std::invalid_argument("Time naming factor must be larger than zero!");
  }
  return factor;
}

/**
 * @brief Gives the checked type of the output times used for the given output
 * type.
 * @param output_type Output type for which the times type should be given.
 * @return times type identifier of the output.
 */
OutputTimesType
OutputReader::ReadOutputTimesType(OutputType const output_type) const {
  return StringToOutputTimesType(DoReadOutputTimesType(output_type));
}

/**
 * @brief Gives the checked interval used to write output for the given ouput
 * type.
 * @param output_type Output type for which the interval should be read.
 * @return interval used for the output.
 */
double OutputReader::ReadOutputInterval(OutputType const output_type) const {
  // Read the time interval and check on consistency
  double const interval(DoReadOutputInterval(output_type));
  if (interval <= 0.0) {
    throw std::invalid_argument("Output interval for " +
                                OutputTypeToString(output_type) +
                                " output must be larger than zero!");
  }
  return interval;
}

/**
 * @brief Gives the checked time stamps in sorted order to be used to write
 * output for the given output type.
 * @param output_type Output type for which the interval should be read.
 * @return positive time stamps in sorted order.
 */
std::vector<double>
OutputReader::ReadOutputTimeStamps(OutputType const output_type) const {
  // Obtain the time stamps and sort them before return
  std::vector<double> time_stamps(DoReadOutputTimeStamps(output_type));
  // If negative elements are present throw error
  if (std::any_of(time_stamps.begin(), time_stamps.end(),
                  [](double const timestamp) { return timestamp < 0.0; })) {
    throw std::invalid_argument("All time stamps for the " +
                                OutputTypeToString(output_type) +
                                " output must be positive or zero!");
  }
  // Sort the time stamps
  std::sort(time_stamps.begin(), time_stamps.end());

  return time_stamps;
}
