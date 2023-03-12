//===------------------------ flux_splitting.h ----------------------------===//
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
#ifndef FLUXSPLITTING_H
#define FLUXSPLITTING_H

#include <iostream>

/**
 * @brief Identifier for flux splitting schemes to be used in Roe Riemann
 * Solver.
 */
enum class FluxSplitting {
  Roe,
  LocalLaxFriedrichs,
  GlobalLaxFriedrichs,
  Roe_M,
  LocalLaxFriedrichs_M
};

/**
 * @brief Provides a string representation of the applied FluxSplitting type in
 * Roe Riemann Solver (used for log file).
 * @param flux_type FluxSplitting identifier.
 * @return String representation for the given FluxSplitting type.
 */
inline std::string FluxSplittingToString(FluxSplitting const flux_type) {
  switch (flux_type) {
  case FluxSplitting::Roe:
    return "Roe";
  case FluxSplitting::LocalLaxFriedrichs:
    return "LocalLaxFriedrichs";
  case FluxSplitting::GlobalLaxFriedrichs:
    return "GlobalLaxFriedrichs";
  case FluxSplitting::Roe_M:
    return "Roe-M";
  case FluxSplitting::LocalLaxFriedrichs_M:
    return "LocalLaxFriedrichs-M";
  default:
    return "ERROR: This flux is not (yet) defined!";
  }
}

#endif // FLUXSPLITTING_H
