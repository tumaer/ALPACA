//===--------------- capillary_pressure_calculator.h ----------------------===//
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
#ifndef CAPILLARY_PRESSURE_CALCULATOR_H
#define CAPILLARY_PRESSURE_CALCULATOR_H

#include "materials/material_manager.h"
#include "topology/node.h"

/**
 * @brief Calculates the capillary pressure.
 */
class CapillaryPressureCalculator {

private:
  /**
   * @brief The surface tension coefficient. Based on it, the capillary pressure
   * is calculated.
   */
  double const surface_tension_coefficient_;

public:
  CapillaryPressureCalculator() = delete;
  explicit CapillaryPressureCalculator(MaterialManager const &material_manager);
  ~CapillaryPressureCalculator() = default;
  CapillaryPressureCalculator(CapillaryPressureCalculator const &) = delete;
  CapillaryPressureCalculator &
  operator=(CapillaryPressureCalculator const &) = delete;
  CapillaryPressureCalculator(CapillaryPressureCalculator &&) = delete;
  CapillaryPressureCalculator &
  operator=(CapillaryPressureCalculator &&) = delete;

  void ComputePressureDifference(
      Node const &node,
      double (&pressure_difference)[CC::TCX()][CC::TCY()][CC::TCZ()]) const;
};

#endif // CAPILLARY_PRESSURE_CALCULATOR_H
