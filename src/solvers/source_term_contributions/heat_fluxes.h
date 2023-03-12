//===--------------------------- heat_fluxes.h ----------------------------===//
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
#ifndef HEAT_FLUXES_H
#define HEAT_FLUXES_H

#include "block_definitions/block.h"
#include "materials/material_manager.h"

/**
 * @brief This class calculates the heat fluxes and adds them to a buffer.
 */
class HeatFluxes {

private:
  MaterialManager const &material_manager_;

public:
  HeatFluxes() = delete;
  explicit HeatFluxes(MaterialManager const &material_manager);
  ~HeatFluxes() = default;
  HeatFluxes(HeatFluxes const &) = delete;
  HeatFluxes(HeatFluxes &&) = delete;
  HeatFluxes &operator=(HeatFluxes const &) = delete;
  HeatFluxes &operator=(HeatFluxes &&) = delete;

  void ComputeFluxes(std::pair<MaterialName const, Block> const &mat_block,
                     double (&heat_fluxes_x)[MF::ANOE()][CC::ICX() + 1]
                                            [CC::ICY() + 1][CC::ICZ() + 1],
                     double (&heat_fluxes_y)[MF::ANOE()][CC::ICX() + 1]
                                            [CC::ICY() + 1][CC::ICZ() + 1],
                     double (&heat_fluxes_z)[MF::ANOE()][CC::ICX() + 1]
                                            [CC::ICY() + 1][CC::ICZ() + 1],
                     double const cell_size) const;
};

#endif // HEAT_FLUXES_H
