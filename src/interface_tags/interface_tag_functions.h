//===------------------- interface_tag_functions.h ------------------------===//
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
#ifndef INTERFACE_TAG_FUNCTIONS_H
#define INTERFACE_TAG_FUNCTIONS_H

#include "user_specifications/compile_time_constants.h"

namespace InterfaceTagFunctions {

void InitializeInternalInterfaceTags(
    std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]);
void SetInternalCutCellTagsFromLevelset(
    double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
    std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]);
void SetTotalInterfaceTagsFromCutCells(
    std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]);
bool TotalInterfaceTagsAreUniform(
    std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]);

} // namespace InterfaceTagFunctions

#endif // INTERFACE_TAG_FUNCTIONS_H
