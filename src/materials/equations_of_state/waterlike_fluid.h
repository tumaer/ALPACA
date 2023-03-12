//===------------------------ waterlike_fluid.h ---------------------------===//
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
#ifndef WATERLIKE_FLUID_H
#define WATERLIKE_FLUID_H

#include "materials/equation_of_state.h"
#include "unit_handler.h"
#include <string>
#include <unordered_map>

/**
 * @brief The WaterlikeFluid class implements a generic waterlike fluid ( aka
 * Tait'S EOS ), i.e. material parameters must be set via caller input.
 * Implementation according to \cite Fedkiw1999a ( some formulas rearranged ).
 */
class WaterlikeFluid : public EquationOfState {
  // variables required for the computation
  double const gamma_;
  double const A_;
  double const B_;
  double const rho0_;

  // functions required from base class
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
  double ComputeSpeedOfSound(double const density,
                             double const pressure) const override;

public:
  WaterlikeFluid() = delete;
  explicit WaterlikeFluid(
      std::unordered_map<std::string, double> const &dimensional_eos_data,
      UnitHandler const &unit_handler);
  virtual ~WaterlikeFluid() = default;
  WaterlikeFluid(WaterlikeFluid const &) = delete;
  WaterlikeFluid &operator=(WaterlikeFluid const &) = delete;
  WaterlikeFluid(WaterlikeFluid &&) = delete;
  WaterlikeFluid &operator=(WaterlikeFluid &&) = delete;

  // function for logging
  std::string GetLogData(unsigned int const indent,
                         UnitHandler const &unit_handler) const;
};

#endif // WATERLIKE_FLUID_H
