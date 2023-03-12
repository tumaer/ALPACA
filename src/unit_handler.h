//===------------------------- unit_handler.h -----------------------------===//
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
#ifndef UNIT_HANDLER_H
#define UNIT_HANDLER_H

#include "enums/unit_type.h"
#include "input_output/log_writer/log_writer.h"
#include <array>
#include <cmath>
#include <vector>

/**
 * @brief The UnitHandler class takes care of the conversions between unitless
 * and non-unitless representations of quantities. Within the Simulation Kernel
 * only non-dimensional quantities are used. In the user input and output,
 * however, values are given in unit representation.
 */
class UnitHandler {
  // Input file reference values and SI units
  double const density_reference_;     // input file
  double const velocity_reference_;    // input file
  double const length_reference_;      // inputfile + SI unit [m]
  double const temperature_reference_; // input file + SI unit [K]
  double const time_reference_;        // SI unit [s]
  double const mass_reference_;        // SI unit [kg]
  // derived reference values for each exisitng unit type
  double const momentum_reference_;
  double const pressure_reference_;
  double const energy_reference_;
  double const viscosity_reference_;
  double const thermal_conductivity_reference_;
  double const surface_tension_coefficient_reference_;

public:
  UnitHandler() = delete;
  explicit UnitHandler(double const density_reference,
                       double const velocity_reference,
                       double const length_reference,
                       double const temperature_reference);
  ~UnitHandler() = default;
  UnitHandler(UnitHandler const &) = delete;
  UnitHandler &operator=(UnitHandler const &) = delete;
  UnitHandler(UnitHandler &&) = delete;
  UnitHandler &operator=(UnitHandler &&) = delete;

  // Dimensionalization and non-dimensionalization functions for a single unit
  // type
  double NonDimensionalizeValue(double const value,
                                UnitType const unit_type) const;
  double DimensionalizeValue(double const value,
                             UnitType const unit_type) const;
  // Dimensionalization and non-dimensionalization for a set of unit types in
  // nominator and denominator
  double
  NonDimensionalizeValue(double const value,
                         std::vector<UnitType> const &nominator_units,
                         std::vector<UnitType> const &denominator_units) const;
  double
  DimensionalizeValue(double const value,
                      std::vector<UnitType> const &nominator_units,
                      std::vector<UnitType> const &denominator_units) const;
};

#endif // UNIT_HANDLER_H
