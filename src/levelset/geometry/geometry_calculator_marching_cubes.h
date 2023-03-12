//===-------------- geomtry_calculator_marching_cubes.h -------------------===//
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
#ifndef GEOMETRY_CALCULATOR_MARCHING_CUBES_H
#define GEOMETRY_CALCULATOR_MARCHING_CUBES_H

#include "geometry_calculator.h"

/**
 * @brief The GeometryCalculatorMarchingCubes class calculates volume fractions
 * and cell face apertures based on the marching-cube algorithm as presented in
 * \cite Lauer2012.
 * @note This class is derived and inherits from the abstract class
 * GeometryCalculator.
 */
class GeometryCalculatorMarchingCubes
    : public GeometryCalculator<GeometryCalculatorMarchingCubes> {

private:
  /**
   * Bool that indicates whether geometric quantities are calculated cell based.
   * The defaults setting is false, and geometric calculations are performed
   * sub-cell based. I.e., A 3-D cell is split into 8 sub cells.
   */
  static constexpr bool cell_based_geometry_calculations_ = false;

public:
  explicit GeometryCalculatorMarchingCubes() = default;
  ~GeometryCalculatorMarchingCubes() = default;
  GeometryCalculatorMarchingCubes(GeometryCalculatorMarchingCubes const &) =
      delete;
  GeometryCalculatorMarchingCubes &
  operator=(GeometryCalculatorMarchingCubes const &) = delete;
  GeometryCalculatorMarchingCubes(GeometryCalculatorMarchingCubes &&) = delete;
  GeometryCalculatorMarchingCubes &
  operator=(GeometryCalculatorMarchingCubes &&) = delete;

  std::array<double, 6> ComputeCellFaceApertureImplementation(
      double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
      unsigned int const i, unsigned int const j, unsigned int const k,
      std::int8_t const material_sign = 1) const;

  double ComputeVolumeFractionImplementation(
      double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
      unsigned int const i, unsigned int const j, unsigned int const k,
      std::int8_t const material_sign = 1) const;
};

#endif // GEOMETRY_CALCULATOR_MARCHING_CUBES_H
