//===---------------------- geometry_calculator.h -------------------------===//
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
#ifndef GEOMETRY_CALCULATOR_H
#define GEOMETRY_CALCULATOR_H

#include "topology/node.h"
#include "user_specifications/compile_time_constants.h"
#include <enums/geometry_settings.h>
#include <vector>

/**
 * @brief The GeometryCalculator class provides functions for basic geometric
 * calculations (normals, curvature, etc.) based on the level-set field and
 * serves as abstract class for more sophisticated geometrical functions
 * (apertures, volume fraction, etc.).
 */
template <typename DerivedGeometryCalculator> class GeometryCalculator {

public:
  GeometryCalculator() = default;
  ~GeometryCalculator() = default;
  GeometryCalculator(GeometryCalculator const &) = delete;
  GeometryCalculator &operator=(GeometryCalculator const &) = delete;
  GeometryCalculator(GeometryCalculator &&) = delete;
  GeometryCalculator &operator=(GeometryCalculator &&) = delete;

  std::array<double, 6> ComputeCellFaceAperture(
      double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
      unsigned int const i, unsigned int const j, unsigned int const k,
      std::int8_t const material_sign = 1) const {
    return static_cast<const DerivedGeometryCalculator *>(this)
        ->ComputeCellFaceApertureImplementation(levelset, i, j, k,
                                                material_sign);
  }

  double ComputeVolumeFraction(
      double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
      unsigned int const i, unsigned int const j, unsigned int const k,
      std::int8_t const material_sign = 1) const {
    return static_cast<const DerivedGeometryCalculator *>(this)
        ->ComputeVolumeFractionImplementation(levelset, i, j, k, material_sign);
  }
};

// box size for subcell and cut-cell computation
static constexpr unsigned int subcell_box_size_x = 3;
static constexpr unsigned int subcell_box_size_y =
    CC::DIM() != Dimension::One ? 3 : 1;
static constexpr unsigned int subcell_box_size_z =
    CC::DIM() == Dimension::Three ? 3 : 1;

static constexpr unsigned int cell_box_size_x = 2;
static constexpr unsigned int cell_box_size_y =
    CC::DIM() != Dimension::One ? 2 : 1;
static constexpr unsigned int cell_box_size_z =
    CC::DIM() == Dimension::Three ? 2 : 1;

void GetLevelsetAtCellCorners(
    double (&cell_corner_levelset)[cell_box_size_x][cell_box_size_y]
                                  [cell_box_size_z],
    double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
    unsigned int const i, unsigned int const j, unsigned int const k);

/**
 * @brief Implements a cut cell criteria.
 * @param levelset The levelset field.
 * @param i The index in x-direction.
 * @param j The index in y-direction.
 * @param k The index in z-direction.
 * @return A bool indicating, whether the given index describes a cut cell.
 * @tparam Different cut cell criteria can b used. Template specifications
 * implement the different criteria.
 */
template <CutCellCriteria>
bool IsCutCell(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
               unsigned int const i, unsigned int const j,
               unsigned int const k);

void ComputeInterfaceCurvature(
    Node const &node, double (&curvature)[CC::TCX()][CC::TCY()][CC::TCZ()]);

void GetLevelsetAtSubcellCorners(
    double (&subcell_corner_levelset)[subcell_box_size_x][subcell_box_size_y]
                                     [subcell_box_size_z],
    double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
    unsigned int const i, unsigned int const j, unsigned int const k);

std::array<double, 3>
GetNormal(double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
          unsigned int const i, unsigned int const j, unsigned int const k,
          std::int8_t const material_sign = 1);

#endif // GEOMETRY_CALCULATOR_H
