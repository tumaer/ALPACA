//===---------------- space_filling_curve_settings.h ----------------------===//
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
#ifndef SPACE_FILLING_CURVE_SETTINGS_H
#define SPACE_FILLING_CURVE_SETTINGS_H

#include "topology/space_filling_curve_index.h"
#include <string>

namespace SpaceFillingCurveSettings {
#if DIMENSION != 3
constexpr auto SfcIndex = LebesgueIndex;
#else
constexpr auto SfcIndex = HilbertIndex;
#endif

/**
 * @brief Gives a string representation of the active space-filling curve.
 */
inline std::string SpaceFillingCurveSelectionString() {
  if (SfcIndex == LebesgueIndex) {
    return "Lebesgue-Curve";
  } else {
    return "Hilbert-Curve";
  }
}
} // namespace SpaceFillingCurveSettings

#endif // SPACE_FILLING_CURVE_SETTINGS_H
