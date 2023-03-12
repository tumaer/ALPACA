//===---------- instantiation_internal_halo_manager.cpp -------------------===//
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
#include "instantiation/halo_manager/instantiation_internal_halo_manager.h"

namespace Instantiation {

/**
 * @brief Instantiates the complete internal halo manager class with the given
 * input classes.
 * @param topology_manager Class providing global (on all ranks) node
 * information.
 * @param tree Tree class providing local (on current rank) node information.
 * @param communication_manager Calls providing communication handling between
 * different ranks
 * @param material_manager material_manager Instance providing initialized
 * material data.
 * @return The fully instantiated InternalHaloManager class.
 */
InternalHaloManager
InstantiateInternalHaloManager(TopologyManager &topology_manager, Tree &tree,
                               CommunicationManager &communication_manager,
                               MaterialManager const &material_manager) {

  return InternalHaloManager(tree, topology_manager, communication_manager,
                             material_manager.GetNumberOfMaterials());
}
} // namespace Instantiation
