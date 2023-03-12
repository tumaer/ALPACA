//===---------------------- axisymmetric_fluxes.h -------------------------===//
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
#ifndef AXISYMMETRIC_FLUXES_H
#define AXISYMMETRIC_FLUXES_H

#include "block_definitions/block.h"

/**
 * This class calculates the additional terms for axisymmetric simulations and
 * adds them to a buffer.
 */
class AxisymmetricFluxes {

public:
  explicit AxisymmetricFluxes() = default;
  ~AxisymmetricFluxes() = default;
  AxisymmetricFluxes(AxisymmetricFluxes const &) = delete;
  AxisymmetricFluxes(AxisymmetricFluxes &&) = delete;
  AxisymmetricFluxes &operator=(AxisymmetricFluxes const &) = delete;
  AxisymmetricFluxes &operator=(AxisymmetricFluxes &&) = delete;

  void ComputeAxisymmetricContributions(
      Block const &block,
      double (&volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()],
      double const cell_size, double const node_origin_x) const;
};

#endif // AXISYMMETRIC_FLUXES_H
