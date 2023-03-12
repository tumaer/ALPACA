//===------------------- geometry_calculator_setup.h ----------------------===//
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
#ifndef GEOMETRY_CALCULATOR_SETUP_H
#define GEOMETRY_CALCULATOR_SETUP_H

#include "geometry_calculator_marching_cubes.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a GeometryCalculator type based on a specified
 * constexpr.
 */
namespace GeometryCalculatorSetup {

/**
 * @brief Function returning the typedef of a GeometryCalculator based on a
 * constexpr template.
 *
 * @tparam GeometryCalculators The constexpr template parameter to specify the
 * exact GeometryCalculator type.
 */
template <GeometryCalculators> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<GeometryCalculators::MarchingCubes> {
  typedef GeometryCalculatorMarchingCubes type;
};

} // namespace GeometryCalculatorSetup

#endif // GEOMETRY_CALCULATOR_SETUP_H
