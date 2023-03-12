//===------------------- multi_phase_manager_setup.h ----------------------===//
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
#ifndef MULTI_PHASE_MANAGER_SETUP_H
#define MULTI_PHASE_MANAGER_SETUP_H

#include "numerous_phase_manager.h"
#include "two_phase_manager.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a MultiPhaseManager type based on a specified
 * constexpr.
 */
namespace MultiPhaseManagerSetup {

/**
 * @brief Function returning the typedef of a MultiPhaseManager based on a
 * constexpr template.
 *
 * @tparam PhaseManagers The constexpr template parameter to specify the exact
 * MultiPhaseManager type.
 */
template <PhaseManagers> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<PhaseManagers::TwoPhase> {
  typedef TwoPhaseManager type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<PhaseManagers::NumerousPhase> {
  typedef NumerousPhaseManager type;
};

} // namespace MultiPhaseManagerSetup

#endif // MULTI_PHASE_MANAGER_SETUP_H
