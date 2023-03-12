//===---------------- parametric_levelset_initializer.h -------------------===//
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
#ifndef PARAMETRIC_LEVELSET_INITIALIZER_H
#define PARAMETRIC_LEVELSET_INITIALIZER_H

#include "initial_condition/levelset_initializer.h"
#include "initial_condition/parametric_variable.h"

/**
 * @brief The ParametricLevelsetInitializer class allows for a levelset
 * initialization based on parametric functions. It contains the levelset
 * computation, the loop structure is in the base class.
 * @note Uses the C++ Mathematical Expression Toolkit Library by Arash Partow,
 * see respective files for License and Copyright information.
 * @note For the parametric levelset some things need to be considered:
 *       1. The number of points for each parametric variable must be large
 * enough to allow proper sampling compared to the cell size.
 *       2. The reference point can be anywhere in the negative material region
 * (for open interfaces, such as quasi 1D in 3D, the reference point should be
 *          placed outside of the domain to ensure that no cells are behind the
 * reference point).
 *       3. The parametric function must be also valid in a small region outside
 * of the domain (halo cells).
 *       4. The parametric function does not contain any undercuts.
 */
class ParametricLevelsetInitializer : public LevelsetInitializer {

  // Data obtained through the constructor
  std::string const parametric_levelset_expession_;
  std::array<ParametricVariable, 2> const parametric_variables_;
  std::array<double, 3> const ref_point_of_negative_levelset_;
  std::vector<std::array<double, 3>> interface_coordinates_;

  // Functions to compute the levelset values
  std::pair<std::array<double, 3>, double>
  FindInterfacePointWithMinimumDistance(
      std::array<double, 3> const &point) const;
  int GetLevelsetSign(std::array<double, 3> const &point,
                      std::array<double, 3> const &interface_point) const;

  // Functions required from base class
  double
  ComputeSignedLevelsetValue(std::array<double, 3> const &point) override;
  std::string GetTypeLogData(unsigned int const indent) const override;

public:
  ParametricLevelsetInitializer() = delete;
  explicit ParametricLevelsetInitializer(
      std::string const &parameteric_expression_string,
      std::array<ParametricVariable, 2> const &parametric_variables,
      std::array<double, 3> const &ref_point_of_positive_levelset,
      std::vector<std::array<double, 6>> const &bounding_boxes,
      std::vector<MaterialName> const &material_names,
      double const node_size_on_level_zero, unsigned int const maximum_level);
  virtual ~ParametricLevelsetInitializer() = default;
  ParametricLevelsetInitializer(ParametricLevelsetInitializer const &) = delete;
  ParametricLevelsetInitializer &
  operator=(ParametricLevelsetInitializer const &) = delete;
  ParametricLevelsetInitializer(ParametricLevelsetInitializer &&) = delete;
  ParametricLevelsetInitializer &
  operator=(ParametricLevelsetInitializer &&) = delete;
};

// Factory functions (outside for testing)
std::vector<std::array<double, 3>> ComputeInterfaceCoordinates(
    std::string const &parametric_expression,
    std::vector<std::string> const &spatial_variable_names,
    std::array<ParametricVariable, 2> const &parametric_variables);

#endif // PARAMETRIC_LEVELSET_INITIALIZER_H
