//===----------------- space_filling_curve_index.h ------------------------===//
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
#ifndef SPACE_FILLING_CURVE_INDEX
#define SPACE_FILLING_CURVE_INDEX

#include "topology/node_id_type.h"
#include <cstdint>

// Space-filling curve index type
using sfcidx_t = std::uint64_t;

sfcidx_t HilbertIndex(nid_t const node_id);

sfcidx_t LebesgueIndex(nid_t const node_id);

#endif // SPACE_FILLING_CURVE_INDEX
