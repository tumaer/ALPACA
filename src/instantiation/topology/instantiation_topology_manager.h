//===--------------- instantiation_topology_manager.h ---------------------===//
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
#ifndef INITIALIZATION_TOPOLOGY_MANAGER_H
#define INITIALIZATION_TOPOLOGY_MANAGER_H

#include "enums/direction_definition.h"
#include "input_output/input_reader.h"
#include "materials/material_manager.h"
#include "topology/topology_manager.h"

/**
 * @brief Defines all instantiation functions required for the topology manager.
 */
namespace Instantiation {

// Factory function to compute the number of blocks per dimension of level zero
std::array<unsigned int, 3> GetNumberOfNodesOnLevelZero(
    MultiResolutionReader const &multi_resolution_reader);

// Functions to return the directions on which periodicity is used
bool IsMaterialPeriodic(
    BoundaryConditionReader const &boundary_condition_reader,
    Direction const direction);
bool IsLevelsetPeriodic(
    BoundaryConditionReader const &boundary_condition_reader,
    Direction const direction);
unsigned int GetActivePeriodicDirections(
    BoundaryConditionReader const &boundary_condition_reader);

// Instantiation function for the topology manager class
TopologyManager
InstantiateTopologyManager(InputReader const &input_reader,
                           MaterialManager const &material_manager);
} // namespace Instantiation

#endif // INITIALIZATION_TOPOLOGY_MANAGER_H
