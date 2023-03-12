//===--------------- functional_levelset_initializer.h --------------------===//
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
#ifndef FUNCTIONAL_LEVELSET_INITIALIZER_H
#define FUNCTIONAL_LEVELSET_INITIALIZER_H

#include "block_definitions/field_interface_definitions.h"
#include "initial_condition/levelset_initializer.h"
#include "user_expression.h"

/**
 * @brief The FunctionalLevelsetInitializer class allows for a levelset
 * initialization based on a levelset function. It contains the levelset
 * computation, the loop structure is in the base class.
 * @note Uses the C++ Mathematical Expression Toolkit Library by Arash Partow,
 * see respective files for License and Copyright information.
 * @note For the functional levelset initializer some things need to be
 * considered:
 *       1. The function must be a signed distance function.
 *       2. The signed distance must be in the magnitude O(1) and must be in the
 * range [-8, 8] in the first cells near the interface.
 */
class FunctionalLevelsetInitializer : public LevelsetInitializer {
  // Member variables for this class only
  std::string const levelset_variable_name_ =
      std::string(IF::InputName(InterfaceDescription::Levelset));
  std::string const levelset_expression_string_;
  std::vector<double> expression_point_ = {0.0, 0.0, 0.0};
  UserExpression const levelset_expression_;

  // Functions required from base class
  double
  ComputeSignedLevelsetValue(std::array<double, 3> const &point) override;
  std::string GetTypeLogData(unsigned int const indent) const override;

public:
  FunctionalLevelsetInitializer() = delete;
  explicit FunctionalLevelsetInitializer(
      std::string const &levelset_expression_string,
      std::vector<std::array<double, 6>> const &bounding_boxes,
      std::vector<MaterialName> const &material_names,
      double const node_size_on_level_zero, unsigned int const maximum_level);
  virtual ~FunctionalLevelsetInitializer() = default;
  FunctionalLevelsetInitializer(FunctionalLevelsetInitializer const &) = delete;
  FunctionalLevelsetInitializer &
  operator=(FunctionalLevelsetInitializer const &) = delete;
  FunctionalLevelsetInitializer(FunctionalLevelsetInitializer &&) = delete;
  FunctionalLevelsetInitializer &
  operator=(FunctionalLevelsetInitializer &&) = delete;
};

#endif // FUNCTIONAL_LEVELSET_INITIALIZER_H
