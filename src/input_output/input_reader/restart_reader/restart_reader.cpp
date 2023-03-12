//===------------------------ restart_reader.cpp --------------------------===//
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
#include "input_output/input_reader/restart_reader/restart_reader.h"

#include <algorithm>
#include <stdexcept>

/**
 * @brief Gives the checked mode used to restore a simulation.
 * @return restore mode identifier of the simulation.
 */
RestoreMode RestartReader::ReadRestoreMode() const {
  return StringToRestoreMode(DoReadRestoreMode());
}

/**
 * @brief Gives the checked filename to be used for the output.
 * @return name of the restore file.
 * @note Filename without any white spaces.
 */
std::string RestartReader::ReadRestoreFilename() const {
  return StringOperations::RemoveSpaces(DoReadRestoreFilename());
}

/**
 * @brief Gives the checked times type used for writing the restart snapshots.
 * @return snapshot times type identifier.
 */
SnapshotTimesType RestartReader::ReadSnapshotTimesType() const {
  return StringToSnapshotTimesType(DoReadSnapshotTimesType());
}

/**
 * @brief Gives the checked interval to be used for writign restart files (in
 * wall seconds).
 * @return interval to be used.
 */
unsigned int RestartReader::ReadSnapshotInterval() const {
  // read and make consistency check
  int const interval(DoReadSnapshotInterval());
  if (interval < 0) {
    throw std::invalid_argument(
        "Snapshot interval for restart files must NOT be below zero!");
  }

  return static_cast<unsigned int>(interval);
}

/**
 * @brief Gives the checked number of snapshots to be kept (for the interval
 * based writing).
 * @return number of snapshots.
 */
unsigned int RestartReader::ReadSnapshotIntervalsToKeep() const {
  // read and make consistency check
  int const keep(DoReadSnapshotIntervalsToKeep());
  if (keep < 0) {
    throw std::invalid_argument(
        "Snapshot intervals to keep for restart files must NOT be below zero!");
  }

  return static_cast<unsigned int>(keep);
}

/**
 * @brief Gives the checked time stamps to be used to write restart snapshots
 * @return timestamps to be used
 */
std::vector<double> RestartReader::ReadSnapshotTimeStamps() const {
  // Obtain the time stamps and sort them before return
  std::vector<double> time_stamps(DoReadSnapshotTimeStamps());
  // If negative elements are present throw error
  if (std::any_of(time_stamps.begin(), time_stamps.end(),
                  [](double const timestamp) { return timestamp < 0.0; })) {
    throw std::invalid_argument(
        "All time stamps for the restart snapshots must be positive or zero!");
  }
  // Sort the time stamps
  std::sort(time_stamps.begin(), time_stamps.end());

  return time_stamps;
}
