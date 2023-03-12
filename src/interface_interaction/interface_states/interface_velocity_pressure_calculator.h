//===------------ interface_velocity_pressure_calculator.h ----------------===//
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
#ifndef INTERFACE_VELOCITY_PRESSURE_CALCULATOR_H
#define INTERFACE_VELOCITY_PRESSURE_CALCULATOR_H

#include "interface_interaction/interface_riemann_solver/interface_riemann_solver_setup.h"
#include "topology/node.h"
#include "user_specifications/numerical_setup.h"

using InterfaceRiemannSolverConcretization =
    InterfaceRiemannSolverSetup::Concretize<interface_riemann_solver>::type;

/**
 * @brief This class calculates the interface pressure and velocity and saves it
 * in the respective buffers of the interface block.
 */
class InterfaceVelocityPressureCalculator {

private:
  InterfaceRiemannSolverConcretization const interface_riemann_solver_;

public:
  InterfaceVelocityPressureCalculator() = delete;
  explicit InterfaceVelocityPressureCalculator(
      MaterialManager const &material_manager);
  ~InterfaceVelocityPressureCalculator() = default;
  InterfaceVelocityPressureCalculator(
      InterfaceVelocityPressureCalculator const &) = delete;
  InterfaceVelocityPressureCalculator &
  operator=(InterfaceVelocityPressureCalculator const &) = delete;
  InterfaceVelocityPressureCalculator(InterfaceVelocityPressureCalculator &&) =
      delete;
  InterfaceVelocityPressureCalculator &
  operator=(InterfaceVelocityPressureCalculator &&) = delete;

  void FillInterfaceVelocityAndPressureBuffer(
      Node &node,
      double (&pressure_difference)[CC::TCX()][CC::TCY()][CC::TCZ()]) const;
};

#endif // INTERFACE_VELOCITY_PRESSURE_CALCULATOR_H
