//===------------------- mesh_generator_utilities.h -----------------------===//
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
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Utility functions for generating the mesh of the current simulation.
 */
namespace MeshGeneratorUtilities {

/**
 * @brief Gives the total number of internal vertices of one single block.
 * @return Number of vertices.
 */
constexpr unsigned int NumberOfInternalVerticesPerBlock() {
  return (CC::ICX() + 1) * (CC::ICY() + 1) * (CC::ICZ() + 1);
}

/**
 * @brief Gives the total number of cells of one single block.
 * @return Number of vertices.
 */
constexpr unsigned int NumberOfTotalVerticesPerBlock() {
  return (CC::TCX() + 1) * (CC::TCY() + 1) * (CC::TCZ() + 1);
}

/**
 * @brief Gives the total number of internal cells of one single block.
 * @return Number of cells.
 */
constexpr unsigned int NumberOfInternalCellsPerBlock() {
  return CC::ICX() * CC::ICY() * CC::ICZ();
}

/**
 * @brief Gives the total number of cells of one single block.
 * @return Number of cells.
 */
constexpr unsigned int NumberOfTotalCellsPerBlock() {
  return CC::TCX() * CC::TCY() * CC::TCZ();
}

/**
 * @brief Gives the cell size for a given block size.
 * @param block_size Size of the block.
 * @return cell_size of the block.
 *
 * @note since only cubic blocks are used, only one direction is required to
 * determine the size. ICX is used since this value is always filled.
 */
constexpr double CellSizeForBlockSize(double const block_size) {
  return block_size / double(CC::ICX());
}
} // namespace MeshGeneratorUtilities
