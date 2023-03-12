//===-------------------- interface_extender_setup.h ----------------------===//
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
#ifndef INTERFACE_EXTENDER_SETUP_H
#define INTERFACE_EXTENDER_SETUP_H

#include "two_phase_interface_extender.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a InterfaceExtender type based on a specified
 * constexpr.
 */
namespace InterfaceExtenderSetup {

/**
 * @brief Function returning the typedef of a InterfaceExtender based on a
 * constexpr template. For each InterfaceFieldType an extension exists.
 *
 * @tparam InterfaceExtenders The constexpr template parameter to specify the
 * exact InterfaceExtender type.
 */
template <InterfaceExtenders> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<InterfaceExtenders::TwoPhase> {
  typedef TwoPhaseInterfaceExtender<InterfaceFieldType::States> type_states;
  typedef TwoPhaseInterfaceExtender<InterfaceFieldType::Parameters>
      type_parameters;
};

} // namespace InterfaceExtenderSetup

#endif // INTERFACE_EXTENDER_SETUP_H
