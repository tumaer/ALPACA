//===---------------------- restart_definitions.h -------------------------===//
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
#ifndef RESTART_DEFINITIONS_H
#define RESTART_DEFINITIONS_H

#include <stdexcept>

#include "utilities/string_operations.h"

/**
 * @brief Identifier for the restore mode used at startup.
 */
enum class RestoreMode { Off, Soft, Forced };

/**
 * @brief The OutputTimesType enum defines the type of the output times to be
 * used for the output. (Off: No output is written) (Interval: A file is written
 * for the given intervals) (Stamps: A file is written for time stamps)
 *        (IntervalStamps: A file is written for interval and time stamps)
 */
enum class SnapshotTimesType { Off, Interval, Stamps, IntervalStamps };

/**
 * @brief Gives the proper OutputWriter type for a given string.
 * @param file_type String that should be converted.
 * @return Output writer type.
 */
inline RestoreMode StringToRestoreMode(std::string const &type) {
  // transform string to upper case without spaces
  std::string const type_upper_case(
      StringOperations::ToUpperCaseWithoutSpaces(type));
  // switch statements cannot be used with strings
  if (type_upper_case == "OFF") {
    return RestoreMode::Off;
  } else if (type_upper_case == "SOFT") {
    return RestoreMode::Soft;
  } else if (type_upper_case == "FORCED") {
    return RestoreMode::Forced;
  } else {
    throw std::logic_error("Restore mode '" + type_upper_case + "' not known!");
  }
}

/**
 * @brief Gives the proper OutputTimes type for a given string.
 * @param times_type String that should be converted.
 * @return Output times type.
 */
inline SnapshotTimesType
StringToSnapshotTimesType(std::string const &times_type) {
  // transform string to upper case without spaces
  std::string const type_upper_case(
      StringOperations::ToUpperCaseWithoutSpaces(times_type));
  // switch statements cannot be used with strings
  if (type_upper_case == "OFF") {
    return SnapshotTimesType::Off;
  } else if (type_upper_case == "INTERVAL") {
    return SnapshotTimesType::Interval;
  } else if (type_upper_case == "STAMPS") {
    return SnapshotTimesType::Stamps;
  } else if (type_upper_case == "INTERVALSTAMPS") {
    return SnapshotTimesType::IntervalStamps;
  } else {
    throw std::logic_error("Restart snaposhot times type '" + type_upper_case +
                           "' not known!");
  }
}

#endif // RESTART_DEFINITIONS_H
