//===---------------------- gravitational_force.h -------------------------===//
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
#ifndef GRAVITATIONAL_FORCE_H
#define GRAVITATIONAL_FORCE_H

#include "block_definitions/block.h"

/**
 * This class calculates the gravitational forces and adds them to a buffer.
 */
class GravitationalForce {

private:
  std::array<double, 3> const gravity_;

public:
  GravitationalForce() = delete;
  explicit GravitationalForce(std::array<double, 3> const gravity);
  ~GravitationalForce() = default;
  GravitationalForce(GravitationalForce const &) = delete;
  GravitationalForce(GravitationalForce &&) = delete;
  GravitationalForce &operator=(GravitationalForce const &) = delete;
  GravitationalForce &operator=(GravitationalForce &&) = delete;

  void ComputeForces(Block const &block,
                     double (&gravity_forces)[MF::ANOE()][CC::ICX()][CC::ICY()]
                                             [CC::ICZ()]) const;
};

#endif // GRAVITATIONAL_FORCE_H
