//===------------------ interface_state_calculator.h ----------------------===//
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
#ifndef INTERFACE_STATE_CALCULATOR_H
#define INTERFACE_STATE_CALCULATOR_H

#include "interface_interaction/interface_states/capillary_pressure_calculator.h"
#include "interface_interaction/interface_states/interface_velocity_pressure_calculator.h"
#include "materials/material_manager.h"
#include "topology/node.h"

/**
 * @brief Calculates the relevant quantities at the interface, e.g. the
 * interface pressure of both materials or the interface velocity.
 */
class InterfaceStateCalculator {
private:
  InterfaceVelocityPressureCalculator const
      interface_velocity_pressure_calculator_;
  CapillaryPressureCalculator const capillary_pressure_calculator_;

public:
  InterfaceStateCalculator() = delete;
  explicit InterfaceStateCalculator(MaterialManager const &material_manager);
  ~InterfaceStateCalculator() = default;
  InterfaceStateCalculator(InterfaceStateCalculator const &) = delete;
  InterfaceStateCalculator &
  operator=(InterfaceStateCalculator const &) = delete;
  InterfaceStateCalculator(InterfaceStateCalculator &&) = delete;
  InterfaceStateCalculator &operator=(InterfaceStateCalculator &&) = delete;

  void ObtainInterfaceStates(Node &node) const;
};

#endif // INTERFACE_STATE_CALCULATOR_H
