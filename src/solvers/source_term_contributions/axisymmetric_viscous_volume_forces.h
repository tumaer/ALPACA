//===-------------- axisymmetric_viscous_volume_forces.h ------------------===//
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
#ifndef AXISYMMETRIC_VISCOUS_VOLUME_FORCES_H
#define AXISYMMETRIC_VISCOUS_VOLUME_FORCES_H

#include <vector>

#include "block_definitions/block.h"
#include "materials/material_manager.h"

/**
 * @brief This class calculates the viscous contribution to axisymmetric forces
 * as described in \cite Meng2016b.
 */
class AxisymmetricViscousVolumeForces {

private:
  MaterialManager const &material_manager_;
  static constexpr unsigned int dim_ =
      2; // only sane configuration; enables unit testing

  void ComputeVelocityGradient(
      double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()],
      double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()],
      double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()],
      double const cell_size,
      double (&velocity_gradient)[CC::TCX()][CC::TCY()][dim_][dim_]) const;

public:
  AxisymmetricViscousVolumeForces() = delete;
  explicit AxisymmetricViscousVolumeForces(
      MaterialManager const &material_manager);
  ~AxisymmetricViscousVolumeForces() = default;
  AxisymmetricViscousVolumeForces(AxisymmetricViscousVolumeForces const &) =
      delete;
  AxisymmetricViscousVolumeForces(AxisymmetricViscousVolumeForces &&) = delete;
  AxisymmetricViscousVolumeForces &
  operator=(AxisymmetricViscousVolumeForces const &) = delete;
  AxisymmetricViscousVolumeForces &
  operator=(AxisymmetricViscousVolumeForces &&) = delete;

  void ComputeForces(
      std::pair<MaterialName const, Block> const &mat_block,
      double (&axisymmetric_viscous_volume_forces)[MF::ANOE()][CC::ICX()]
                                                  [CC::ICY()][CC::ICZ()],
      double const cell_size, double const node_origin_x) const;
};

#endif // AXISYMMETRIC_VISCOUS_VOLUME_FORCES_H
