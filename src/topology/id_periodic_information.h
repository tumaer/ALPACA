//===------------------- id_periodic_information.h ------------------------===//
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
#ifndef ID_PERIODIC_INFORMATION_H
#define ID_PERIODIC_INFORMATION_H

#include "boundary_condition/boundary_specifications.h"
#include "topology/node_id_type.h"
#include <array>

/**
 * @brief Unique Identifier for the location of the periodic boundaries.
 * Designed to use as bitmask for finding if the x, y or z axis are periodic.
 */
enum PeriodicBoundariesLocations : unsigned int {
  EastWest = 1 << 0,
  NorthSouth = 1 << 1,
  TopBottom = 1 << 2
};

nid_t GetPeriodicNeighborId(
    nid_t const id, BoundaryLocation const location,
    std::array<unsigned int, 3> const level_zero_blocks_xyz,
    unsigned int const active_periodic_locations);
bool PeriodicIsExternalBoundary(
    BoundaryLocation const location, nid_t const id,
    std::array<unsigned int, 3> const level_zero_blocks_xyz,
    unsigned int const active_periodic_locations);

#endif // ID_PERIODIC_INFORMATION_H
