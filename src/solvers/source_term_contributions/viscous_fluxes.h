//===-------------------------- viscous_fluxes.h --------------------------===//
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
#ifndef VISCOUS_FLUXES_H
#define VISCOUS_FLUXES_H

#include <vector>

#include "block_definitions/block.h"
#include "materials/material_manager.h"

/**
 * @brief This class calculates the viscous source terms and adds them to a flux
 * buffer.
 */
class ViscousFluxes {

private:
  MaterialManager const &material_manager_;

  void ComputeTauFluxes(
      double const (
          &velocity_gradient_at_cell_faces)[CC::ICX() + 1][CC::ICY() + 1]
                                           [CC::ICZ() + 1][DTI(CC::DIM())]
                                           [DTI(CC::DIM())][DTI(CC::DIM())],
      std::vector<double> const viscosity,
      double (&tau)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI(CC::DIM())]
                   [DTI(CC::DIM())]) const;

  void ComputeTauFluxes(
      double const (
          &velocity_gradient_at_cell_faces)[CC::ICX() + 1][CC::ICY() + 1]
                                           [CC::ICZ() + 1][DTI(CC::DIM())]
                                           [DTI(CC::DIM())][DTI(CC::DIM())],
      double const (
          &shear_viscosity_at_cell_faces)[CC::ICX() + 1][CC::ICY() + 1]
                                         [CC::ICZ() + 1][DTI(CC::DIM())],
      double const bulk_viscosity,
      double (&tau)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI(CC::DIM())]
                   [DTI(CC::DIM())]) const;

public:
  ViscousFluxes() = delete;
  explicit ViscousFluxes(MaterialManager const &material_manager);
  ~ViscousFluxes() = default;
  ViscousFluxes(ViscousFluxes const &) = delete;
  ViscousFluxes(ViscousFluxes &&) = delete;
  ViscousFluxes &operator=(ViscousFluxes const &) = delete;
  ViscousFluxes &operator=(ViscousFluxes &&) = delete;

  void ComputeFluxes(std::pair<MaterialName const, Block> const &mat_block,
                     double (&dissipative_flux_x)[MF::ANOE()][CC::ICX() + 1]
                                                 [CC::ICY() + 1][CC::ICZ() + 1],
                     double (&dissipative_flux_y)[MF::ANOE()][CC::ICX() + 1]
                                                 [CC::ICY() + 1][CC::ICZ() + 1],
                     double (&dissipative_flux_z)[MF::ANOE()][CC::ICX() + 1]
                                                 [CC::ICY() + 1][CC::ICZ() + 1],
                     double cell_size) const;
};

#endif // VISCOUS_FLUXES_H
