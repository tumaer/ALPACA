//===--------------------------- unit_type.h ------------------------------===//
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
#ifndef UNIT_TYPE_H
#define UNIT_TYPE_H

/**
 * @brief Identifier for unit types used in (non-)dimensionalization of physical
 * quantities.
 */
enum class UnitType {
  Unitless,
  // SI units
  Length,
  Time,
  Temperature,
  Mass,
  // conservatives/prime states
  Density,
  Velocity,
  Momentum,
  Energy,
  Pressure,
  // material parameter
  Viscosity,
  ThermalConductivity,
  // interface parameter
  SurfaceTensionCoefficient
};

#endif // UNIT_TYPE_H
