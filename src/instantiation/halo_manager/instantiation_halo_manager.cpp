//===---------------- instantiation_halo_manager.cpp ----------------------===//
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
#include "instantiation/halo_manager/instantiation_halo_manager.h"

namespace Instantiation {

/**
 * @brief Instantiates the complete internal halo manager class with the given
 * input classes.
 * @param topology_manager Class providing global (on all ranks) node
 * information.
 * @param tree Tree class providing local (on current rank) node information.
 * @param external_halo_manager Instance for handling external boundary
 * conditions.
 * @param internal_halo_manager Instance for handling internal boundary
 * conditions.
 * @param communication_manager Calls providing communication handling between
 * different ranks.
 * @return The fully instantiated HaloManager class.
 */
HaloManager
InstantiateHaloManager(TopologyManager const &topology_manager, Tree &tree,
                       ExternalHaloManager const &external_halo_manager,
                       InternalHaloManager &internal_halo_manager,
                       CommunicationManager &communication_manager) {

  // return the initializedhalo manager
  return HaloManager(tree, external_halo_manager, internal_halo_manager,
                     communication_manager, topology_manager.GetMaximumLevel());
}
} // namespace Instantiation
