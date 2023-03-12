//===------------------- ghost_fluid_extender_setup.h ---------------------===//
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
#ifndef GHOST_FLUID_EXTENDER_SETUP_H
#define GHOST_FLUID_EXTENDER_SETUP_H

#include "levelset/multi_phase_manager/ghost_fluid_extender/fedkiw_iterative_ghost_fluid_extender.h"
#include "levelset/multi_phase_manager/ghost_fluid_extender/upwind_iterative_ghost_fluid_extender.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a GhostFluidExtender type based on a specified
 * constexpr.
 */
namespace GhostFluidExtenderSetup {

/**
 * @brief Function returning the typedef of a GhostFluidExtender based on a
 * constexpr template. For each MaterialFieldType an extension exists.
 *
 * @tparam GhostFluidExtenders The constexpr template parameter to specify the
 * exact GhostFluidExtender type.
 */
template <Extenders> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<Extenders::Fedkiw> {
  typedef FedkiwGhostFluidExtender<MaterialFieldType::Conservatives>
      type_conservatives;
  typedef FedkiwGhostFluidExtender<MaterialFieldType::PrimeStates>
      type_primestates;
  typedef FedkiwGhostFluidExtender<MaterialFieldType::Parameters>
      type_parameters;
};

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<Extenders::Upwind> {
  typedef UpwindGhostFluidExtender<MaterialFieldType::Conservatives>
      type_conservatives;
  typedef UpwindGhostFluidExtender<MaterialFieldType::PrimeStates>
      type_primestates;
  typedef UpwindGhostFluidExtender<MaterialFieldType::Parameters>
      type_parameters;
};
} // namespace GhostFluidExtenderSetup

#endif // GHOST_FLUID_EXTENDER_SETUP_H
