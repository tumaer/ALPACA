//===----------------- instantiation_halo_manager.h -----------------------===//
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
#ifndef INSTANTIATION_HALO_MANAGER_H
#define INSTANTIATION_HALO_MANAGER_H

#include "halo_manager.h"
#include "topology/topology_manager.h"

/**
 * @brief Defines all instantiation functions required for the halo manager.
 */
namespace Instantiation {

// Instantiation function of the halo manager
HaloManager
InstantiateHaloManager(TopologyManager const &topology_manager, Tree &tree,
                       ExternalHaloManager const &external_halo_manager,
                       InternalHaloManager &internal_halo_manager,
                       CommunicationManager &communication_manager);
} // namespace Instantiation

#endif // INSTANTIATION_HALO_MANAGER_H
