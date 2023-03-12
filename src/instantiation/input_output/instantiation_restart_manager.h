//===----------------- instantiation_restart_manager.h --------------------===//
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
#ifndef INSTANTIATION_RESTART_MANAGER_H
#define INSTANTIATION_RESTART_MANAGER_H

#include <memory>
#include <vector>

#include "input_output/restart_manager.h"

/**
 * @brief Defines all instantiation functions required for the restart manager.
 */
namespace Instantiation {

// Instantiation function for the input_output manager
RestartManager InstantiateRestartManager(TopologyManager &topology_manager,
                                         Tree &tree,
                                         UnitHandler const &unit_handler);
} // namespace Instantiation

#endif // INSTANTIATION_RESTART_MANAGER_H
