//===----------------------- stiffened_gas_safe.h -------------------------===//
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
#ifndef STIFFENED_GAS_SAFE_H
#define STIFFENED_GAS_SAFE_H

#include "materials/equation_of_state.h"
#include "unit_handler.h"
#include <limits>
#include <unordered_map>

/**
 * @brief The StiffenedGasSafe class implements a generic stiffened gas, i.e.
 * material parameters must be set via caller input. Safe version, i.e. uses
 * epsilons to avoid floating-point errors.
 */
class StiffenedGasSafe : public EquationOfState {
  // parameter required for the computataion
  double const epsilon_ = std::numeric_limits<double>::epsilon();
  double const gamma_;
  double const background_pressure_;

  // functions required from the base class
  double ComputePressure(double const mass, double const momentum_x,
                         double const momentum_y, double const momentum_z,
                         double const energy) const override;
  double ComputeEnthalpy(double const mass, double const momentum_x,
                         double const momentum_y, double const momentum_z,
                         double const energy) const override;
  double ComputeEnergy(double const density, double const velocity_x,
                       double const velocity_y, double const velocity_z,
                       double const pressure) const override;
  double GetGruneisen() const override;
  double ComputePsi(double const pressure,
                    double const one_density) const override;
  double GetGamma() const override;
  double GetB() const override;
  double ComputeSpeedOfSound(double const density,
                             double const pressure) const override;

public:
  StiffenedGasSafe() = delete;
  explicit StiffenedGasSafe(
      std::unordered_map<std::string, double> const &dimensional_eos_data,
      UnitHandler const &unit_handler);
  virtual ~StiffenedGasSafe() = default;
  StiffenedGasSafe(StiffenedGasSafe const &) = delete;
  StiffenedGasSafe &operator=(StiffenedGasSafe const &) = delete;
  StiffenedGasSafe(StiffenedGasSafe &&) = delete;
  StiffenedGasSafe &operator=(StiffenedGasSafe &&) = delete;

  // function for logging
  std::string GetLogData(unsigned int const indent,
                         UnitHandler const &unit_handler) const;
};

#endif // STIFFENED_GAS_SAFE_H
