//===--------------------- output_definitions.h ---------------------------===//
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
#ifndef OUTPUT_DEFINITIONS_H
#define OUTPUT_DEFINITIONS_H

#include <stdexcept>
#include <vector>

#include "utilities/string_operations.h"

/**
 * @brief The OutputTimesType enum defines the type of the output times to be
 * used for the output. (Off: No output is written). (Interval: A file is
 * written for the given intervals). (Stamps: A file is written for time
 * stamps). (IntervalStamps: A file is written for the given intervals and
 * period).
 */
enum class OutputTimesType { Off, Interval, Stamps, IntervalStamps };

/**
 * @brief The OutputType defines the type of output which is written (standard,
 * interface and debug).
 */
// AB 2020-03-23 Do not change underlying type and indices. Used for mapping of
// correct array position.
enum class OutputType : unsigned short {
  Standard = 0,
  Interface = 1,
  Debug = 2
};

/**
 * @brief Converts an output type identifier to a (C++11 standard compliant, i.
 * e. positive) array index. "OTTI = Output Type To Index"
 * @param ot The output type identifier
 * @return The index.
 */
constexpr std::underlying_type<OutputType>::type OTTI(OutputType const ot) {
  return static_cast<typename std::underlying_type<OutputType>::type>(ot);
}

/**
 * @brief Converts the OutputType to its corresponding string (for logging).
 * @param type The Output type identifier.
 * @return String to be used.
 */
inline std::string OutputTypeToString(OutputType const type) {

  switch (type) {
  case OutputType::Standard: {
    return "Standard";
  }
  case OutputType::Interface: {
    return "Interface";
  }
  case OutputType::Debug: {
    return "Debug";
  }
  default: {
    throw std::logic_error("Output type is not known!");
  }
  }
}

/**
 * @brief Gives the proper OutputTimes type for a given string.
 * @param times_type String that should be converted.
 * @return Output times type.
 */
inline OutputTimesType StringToOutputTimesType(std::string const &times_type) {
  // transform string to upper case without spaces
  std::string const type_upper_case(
      StringOperations::ToUpperCaseWithoutSpaces(times_type));
  // switch statements cannot be used with strings
  if (type_upper_case == "OFF") {
    return OutputTimesType::Off;
  } else if (type_upper_case == "INTERVAL") {
    return OutputTimesType::Interval;
  } else if (type_upper_case == "STAMPS") {
    return OutputTimesType::Stamps;
  } else if (type_upper_case == "INTERVALSTAMPS" ||
             type_upper_case == "STAMPSINTERVAL") {
    return OutputTimesType::IntervalStamps;
  } else {
    throw std::logic_error("Output times type '" + type_upper_case +
                           "' not known!");
  }
}

#endif // OUTPUT_DEFINITIONS_H
