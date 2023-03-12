//===-------------- instantiation_restart_manager.cpp ---------------------===//
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
#include "instantiation/input_output/instantiation_restart_manager.h"

namespace Instantiation {

/**
 * @brief Instantiates the complete restart manager class with the given input
 * classes.
 * @param topology_manager Class providing global (on all ranks) node
 * information.
 * @param tree Tree class providing local (on current rank) node information.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @return The fully instantiated RestartManager class as pointer (allows
 * movements of it).
 */
RestartManager InstantiateRestartManager(TopologyManager &topology_manager,
                                         Tree &tree,
                                         UnitHandler const &unit_handler) {

  // return the fully initialied restart manager
  return RestartManager(unit_handler, topology_manager, tree,
                        topology_manager.GetMaximumLevel());
}
} // namespace Instantiation
