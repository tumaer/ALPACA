//===--------------------- axisymmetric_fluxes.cpp ------------------------===//
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
#include "axisymmetric_fluxes.h"

/**
 * @brief Computes source terms for axisymmetric simulations. Terms are set
 * according to \cite Adami2016.
 * @param block Block of the considered phase.
 * @param volume_forces Reference to array of volume forces increments to be
 * filled here (indirect return parameter).
 */
void AxisymmetricFluxes::ComputeAxisymmetricContributions(
    Block const &block,
    double (&volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()],
    double const cell_size, double const node_origin_x) const {

  double const(&velocity_x)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      block.GetPrimeStateBuffer(PrimeState::VelocityX);
  double const(&pressure)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      block.GetPrimeStateBuffer(PrimeState::Pressure);
  double const(&momentum_x)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      block.GetAverageBuffer(Equation::MomentumX);
  // direct use of y-momentum buffer allowed since axisymmetric is only used in
  // 2D
  double const(&momentum_y)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      block.GetAverageBuffer(Equation::MomentumY);
  double const(&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      block.GetAverageBuffer(Equation::Energy);

  for (unsigned int i = 0; i < CC::ICX(); ++i) {
    double const one_radius =
        1.0 / (node_origin_x + ((double(i) + 0.5)) * cell_size);
    unsigned int const index_i = i + CC::FICX();
    for (unsigned int j = 0; j < CC::ICY(); ++j) {
      unsigned int const index_j = j + CC::FICY();

      volume_forces[ETI(Equation::Mass)][i][j][0] -=
          momentum_x[index_i][index_j][0] * one_radius;
      if constexpr (MF::IsEquationActive(Equation::Energy)) {
        volume_forces[ETI(Equation::Energy)][i][j][0] -=
            velocity_x[index_i][index_j][0] *
            (energy[index_i][index_j][0] + pressure[index_i][index_j][0]) *
            one_radius;
      }
      volume_forces[ETI(Equation::MomentumX)][i][j][0] -=
          velocity_x[index_i][index_j][0] * momentum_x[index_i][index_j][0] *
          one_radius;
      volume_forces[ETI(Equation::MomentumY)][i][j][0] -=
          velocity_x[index_i][index_j][0] * momentum_y[index_i][index_j][0] *
          one_radius;
    }
  }
}
