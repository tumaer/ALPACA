//===------------ instantiation_internal_halo_manager.h -------------------===//
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
#ifndef INSTANTIATION_INTERNAL_HALO_MANAGER_H
#define INSTANTIATION_INTERNAL_HALO_MANAGER_H

#include "communication/internal_halo_manager.h"
#include "materials/material_manager.h"

/**
 * @brief Defines all instantiation functions required for the internal halo
 * manager.
 */
namespace Instantiation {

// Instantiation function of the internal halo manager
InternalHaloManager
InstantiateInternalHaloManager(TopologyManager &topology_manager, Tree &tree,
                               CommunicationManager &communication_manager,
                               MaterialManager const &material_manager);
} // namespace Instantiation

#endif // INSTANTIATION_INTERNAL_HALO_MANAGER_H
