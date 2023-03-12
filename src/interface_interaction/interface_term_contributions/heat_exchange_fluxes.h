//===--------------------- heat_exchange_fluxes.h -------------------------===//
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
#ifndef HEAT_EXCHANGE_FLUXES_H
#define HEAT_EXCHANGE_FLUXES_H

#include "topology/node.h"

/**
 * @brief This class calculates heat fluxes over the interface and adds them to
 * the right-hand side.
 */
class HeatExchangeFluxes {

private:
  double const thermal_conductivity_positive_;
  double const thermal_conductivity_negative_;

  static constexpr double epsilon_ = std::numeric_limits<double>::epsilon();

  void ComputeRealMaterialTemperature(
      Node const &node,
      double (
          &real_material_temperature)[CC::TCX()][CC::TCY()][CC::TCZ()]) const;

public:
  HeatExchangeFluxes() = delete;
  explicit HeatExchangeFluxes(double const thermal_conductivity_positive,
                              double const thermal_conductivity_negative);
  ~HeatExchangeFluxes() = default;
  HeatExchangeFluxes(HeatExchangeFluxes const &) = delete;
  HeatExchangeFluxes &operator=(HeatExchangeFluxes const &) = delete;
  HeatExchangeFluxes(HeatExchangeFluxes &&) = delete;
  HeatExchangeFluxes &operator=(HeatExchangeFluxes &&) = delete;

  void ComputeInterfaceFluxes(
      Node &node,
      double const (
          &delta_aperture_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3]) const;
};

#endif // HEAT_EXCHANGE_FLUXES_H
