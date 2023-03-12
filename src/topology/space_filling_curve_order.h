//===------------------ space_filling_curve_order.h -----------------------===//
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
#ifndef SPACE_FILLING_CURVE_ORDER_H
#define SPACE_FILLING_CURVE_ORDER_H

#include "topology/node_id_type.h"
#include "topology/space_filling_curve_index.h"
#include "user_specifications/space_filling_curve_settings.h"
#include <algorithm>
#include <iterator>
#include <vector>

/**
 * @brief Sorts a list of node ids according the the provided space-filling
 * curve .
 * @param ids_to_sort The unsorted indices, which get sorted by this function.
 * Indirect return parameter.
 * @param index The index function of the respective space-filling curve.
 */
template <typename SpaceFillingCurveIndexFunction =
              decltype(SpaceFillingCurveSettings::SfcIndex)>
void OrderNodeIdsBySpaceFillingCurve(std::vector<nid_t> &ids_to_sort,
                                     SpaceFillingCurveIndexFunction index =
                                         SpaceFillingCurveSettings::SfcIndex) {
  std::sort(
      std::begin(ids_to_sort), std::end(ids_to_sort),
      [&index](nid_t const a, nid_t const b) { return index(a) < index(b); });
}

#endif // SPACE_FILLING_CURVE_ORDER_H
