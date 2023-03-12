//===-------------------- time_control_reader.cpp -------------------------===//
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
#include "input_output/input_reader/time_control_reader/time_control_reader.h"

#include <stdexcept>

/**
 * @brief Gives the checked start time from the input.
 * @return start time of the simulation.
 */
double TimeControlReader::ReadStartTime() const {
  // Read the start time and check if it is larger than zero
  double const start_time(DoReadStartTime());
  if (start_time < 0.0) {
    throw std::invalid_argument("Start time must be larger than zero!");
  }

  return start_time;
}

/**
 * @brief Gives the checked end time from the input.
 * @return end time of the simulation.
 */
double TimeControlReader::ReadEndTime() const {
  // Read the start and end time and check if end time is larger than start time
  // (no check for negative required since then start time would throw error)
  double const start_time(DoReadStartTime());
  double const end_time(DoReadEndTime());
  if (end_time < start_time) {
    throw std::invalid_argument("End time must be larger than start time!");
  }

  return end_time;
}

/**
 * @brief Gives the checked CFL number from the input.
 * @return CFL umber of the simulation.
 */
double TimeControlReader::ReadCFLNumber() const {
  // Read the start time and check if it is larger than one or below zero
  double const cfl_number(DoReadCFLNumber());
  if (cfl_number < 0.0 || cfl_number > 1.0) {
    throw std::invalid_argument("CFL number must be between 0 and 1 !");
  }

  return cfl_number;
}
