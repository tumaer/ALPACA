//===---------------------- parametric_variable.h -------------------------===//
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
#ifndef PARAMETRIC_VARIABLE_H
#define PARAMETRIC_VARIABLE_H

#include "utilities/string_operations.h"
#include <string>

struct ParametricVariable {
  std::string name = "";
  std::uint64_t points = 1;
  double start = 0.0;
  double end = 0.0;
  double delta = 0.0;

  ParametricVariable() = default;
  ~ParametricVariable() = default;
  ParametricVariable(ParametricVariable const &) = default;
  ParametricVariable &operator=(ParametricVariable const &) = default;
  ParametricVariable(ParametricVariable &&) = default;
  ParametricVariable &operator=(ParametricVariable &&) = default;

  /**
   * @brief Constructs a parametric variable.
   * @param name The name of the variable.
   * @param start The start point of the variable.
   * @param end The end point of the variable.
   * @param points The number of points between start and end (inclusive both).
   */
  explicit ParametricVariable(std::string const var_name,
                              double const var_start, double const var_end,
                              std::uint64_t const var_points)
      : name(var_name), points(var_points > 0 ? var_points : 1),
        start(CreateStart(var_start, var_end)), end(var_end),
        delta(CreateDelta(start, end, points)) {
    // Empty besides initializer list
  }

  /**
   * @brief Gives the start value and makes sanity check that it is larger than
   * the end value.
   * @return The start value.
   */
  double CreateStart(double const start, double const end) {
    if (end < start) {
      throw std::invalid_argument("For the parametric variable the start value "
                                  "must be below the end value.");
    }
    return start;
  }

  /**
   * @brief Gives the delta value for the given start, end and points.
   * @return The delta value.
   */
  double CreateDelta(double const start, double const end,
                     std::uint64_t const points) {
    if (points > 1) {
      double const delta = (end - start) / double(points - 1);
      return delta > 0 ? delta : 0.0;
    } else {
      return 0.0;
    }
  }

  /**
   * @brief Logs the data for this variable.
   * @return The log string.
   */
  std::string GetLogData(unsigned int const indent) const {
    std::string log_string =
        StringOperations::Indent(indent) + "Name: " + name + "\n";
    log_string += StringOperations::Indent(indent + 2) + "Start : " +
                  StringOperations::ToScientificNotationString(start, 2) + "\n";
    log_string += StringOperations::Indent(indent + 2) + "End   : " +
                  StringOperations::ToScientificNotationString(end, 2) + "\n";
    log_string += StringOperations::Indent(indent + 2) + "Delta : " +
                  StringOperations::ToScientificNotationString(delta, 2) + "\n";
    log_string += StringOperations::Indent(indent + 2) +
                  "Points: " + std::to_string(points) + "\n";
    return log_string;
  }
};

#endif // PARAMETRIC_VARIABLE_H
