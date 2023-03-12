//===------------------ instantiation_unit_handler.h ----------------------===//
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
#ifndef INITIALIZATION_UNIT_HANDLER_H
#define INITIALIZATION_UNIT_HANDLER_H

#include "input_output/input_reader.h"
#include "unit_handler.h"

/**
 * @brief Defines all instantiation functions required for the unit handler.
 */
namespace Instantiation {

// Instantiation function of the unit handler
UnitHandler InstantiateUnitHandler(InputReader const &input_reader);
} // namespace Instantiation

#endif // INITIALIZATION_UNIT_HANDLER_H
