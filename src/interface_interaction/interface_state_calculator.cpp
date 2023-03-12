//===---------------- interface_state_calculator.cpp ----------------------===//
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
#include "interface_state_calculator.h"
#include "levelset/multi_phase_manager/two_phase_manager.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief A constructor for the InterfaceStateCalculator.
 * @param material_manager A material manager giving information about the
 * materials present in the simulation.
 */
InterfaceStateCalculator::InterfaceStateCalculator(
    MaterialManager const &material_manager)
    : interface_velocity_pressure_calculator_(material_manager),
      capillary_pressure_calculator_(material_manager) {
  // Empty besides initializer list.
}

/**
 * @brief Computes the interface quantities for a node.
 * @param node The node for which the interface quantities are calculated.
 * @note This function should only be called for multi-phase nodes
 * that have a level-set block. Sanity check whether the input is a multi-phase
 * node is not done!
 */
void InterfaceStateCalculator::ObtainInterfaceStates(Node &node) const {

  double pressure_difference[CC::TCX()][CC::TCY()][CC::TCZ()];
  for (unsigned int i = 0; i < CC::TCX(); ++i) {
    for (unsigned int j = 0; j < CC::TCY(); ++j) {
      for (unsigned int k = 0; k < CC::TCZ(); ++k) {
        pressure_difference[i][j][k] = 0.0;
      }
    }
  }

  if (CC::CapillaryForcesActive() && CC::DIM() != Dimension::One) {
    capillary_pressure_calculator_.ComputePressureDifference(
        node, pressure_difference);
  }
  interface_velocity_pressure_calculator_
      .FillInterfaceVelocityAndPressureBuffer(node, pressure_difference);
}
