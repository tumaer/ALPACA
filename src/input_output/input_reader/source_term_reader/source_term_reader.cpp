//===--------------------- source_term_reader.cpp -------------------------===//
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
#include "input_output/input_reader/source_term_reader/source_term_reader.h"

#include "enums/direction_definition.h"

/**
 * @brief Gives the gravity in a certain direction.
 * @param direction Direction in which the gravity should be read.
 * @return gravity in that direction.
 *
 * @note No check on negative values is required. This would prevent the
 * declaration of negative accelerations.
 */
double SourceTermReader::ReadGravity(Direction const direction) const {
  return DoReadGravity(direction);
}
