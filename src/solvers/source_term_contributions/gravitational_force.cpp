//===---------------------- gravitational_force.cpp -----------------------===//
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
#include "gravitational_force.h"

#include "utilities/index_transformations.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief The default constructor of the class.
 * @param gravity The vector containing the components of the gravitational
 * force.
 */
GravitationalForce::GravitationalForce(std::array<double, 3> const gravity)
    : gravity_(gravity) {
  // Empty besides initializer list
}

/**
 * @brief Computes increments for cell averages due to gravity.
 * @param block Block of the considered phase.
 * @param gravity_forces Reference to array of volume forces increments to be
 * filled here (indirect return parameter).
 */
void GravitationalForce::ComputeForces(
    Block const &block,
    double (
        &gravity_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()]) const {

  Conservatives const &conservatives = block.GetAverageBuffer();
  double const(&density)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      block.GetAverageBuffer(Equation::Mass);

  for (unsigned int i = 0; i < CC::ICX(); ++i) {
    for (unsigned int j = 0; j < CC::ICY(); ++j) {
      for (unsigned int k = 0; k < CC::ICZ(); ++k) {

        std::array<unsigned int, 3> const indices = {BIT::I2TX(i), BIT::I2TY(j),
                                                     BIT::I2TZ(k)};

        // Add up to volume forces
        gravity_forces[ETI(Equation::Mass)][i][j][k] += 0.0;
        if constexpr (MF::IsEquationActive(Equation::Energy)) {
          gravity_forces[ETI(Equation::Energy)][i][j][k] +=
              DimensionAwareConsistencyManagedSum(
                  gravity_[0] * conservatives[Equation::MomentumX][indices[0]]
                                             [indices[1]][indices[2]],
                  CC::DIM() != Dimension::One
                      ? gravity_[1] *
                            conservatives[Equation::MomentumY][indices[0]]
                                         [indices[1]][indices[2]]
                      : 0.0,
                  CC::DIM() == Dimension::Three
                      ? gravity_[2] *
                            conservatives[Equation::MomentumZ][indices[0]]
                                         [indices[1]][indices[2]]
                      : 0.0);
        }
        gravity_forces[ETI(Equation::MomentumX)][i][j][k] +=
            gravity_[0] * density[indices[0]][indices[1]][indices[2]];
        if constexpr (MF::IsEquationActive(Equation::MomentumY)) {
          gravity_forces[ETI(Equation::MomentumY)][i][j][k] +=
              gravity_[1] * density[indices[0]][indices[1]][indices[2]];
        }
        if constexpr (MF::IsEquationActive(Equation::MomentumZ)) {
          gravity_forces[ETI(Equation::MomentumZ)][i][j][k] +=
              gravity_[2] * density[indices[0]][indices[1]][indices[2]];
        }
      }
    }
  }
}
