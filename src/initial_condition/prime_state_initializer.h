//===------------------- prime_state_initializer.h ------------------------===//
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
#ifndef PRIME_STATE_INITIALIZER_H
#define PRIME_STATE_INITIALIZER_H

#include <string>
#include <vector>

#include "block_definitions/field_material_definitions.h"
#include "materials/material_definitions.h"
#include "topology/node_id_type.h"
#include "unit_handler.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief The PrimeStateInitializer class allows for a prime state
 * initialization.
 * @note Uses the C++ Mathematical Expression Toolkit Library by Arash Partow,
 * see respective files for License and Copyright information.
 */
class PrimeStateInitializer {
  // Instance for dimensionalization and non-dimensionalization
  UnitHandler const &unit_handler_;

  // Member variable providing the user expression of input data for materials
  std::vector<std::string> const prime_state_expression_strings_;
  std::vector<std::string> const prime_state_variable_names_;
  std::vector<std::string> const spatial_variable_names_ = {"x", "y", "z"};

  // Additional required variables
  double const dimensionalized_node_size_on_level_zero_;

public:
  PrimeStateInitializer(PrimeStateInitializer const &) = delete;
  explicit PrimeStateInitializer(
      std::vector<std::string> const &prime_state_expression_strings,
      std::vector<std::string> const &prime_state_variable_names,
      double const dimensionalized_node_size_on_level_zero,
      UnitHandler const &unit_handler);
  ~PrimeStateInitializer() = default;
  PrimeStateInitializer &operator=(PrimeStateInitializer const &) = delete;
  PrimeStateInitializer(PrimeStateInitializer &&) = delete;
  PrimeStateInitializer &operator=(PrimeStateInitializer &&) = delete;

  // Public function that can be called from outside
  void GetInitialPrimeStates(
      nid_t const node_id, MaterialName const material,
      double (
          &initial_values)[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()]) const;
};

#endif // PRIME_STATE_INITIALIZER_H
