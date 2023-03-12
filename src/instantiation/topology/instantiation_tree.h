//===----------------------- instantiation_tree.h -------------------------===//
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
#ifndef INITIALIZATION_TREE_H
#define INITIALIZATION_TREE_H

#include "input_output/input_reader.h"
#include "topology/tree.h"
#include "unit_handler.h"

/**
 * @brief Defines all instantiation functions required for the tree.
 */
namespace Instantiation {

// Instantiation function of the tree
Tree InstantiateTree(InputReader const &input_reader,
                     TopologyManager &topology_manager,
                     UnitHandler const &unit_handler);
} // namespace Instantiation

#endif // INITIALIZATION_TREE_H
