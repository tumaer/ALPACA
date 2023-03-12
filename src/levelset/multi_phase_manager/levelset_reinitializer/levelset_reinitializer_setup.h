//===------------------- levelset_reinitializer_setup.h -------------------===//
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
#ifndef LEVELSET_REINITIALIZER_SETUP_H
#define LEVELSET_REINITIALIZER_SETUP_H

#include "min_iterative_levelset_reinitializer.h"
#include "user_specifications/numerical_setup.h"
#include "weno_iterative_levelset_reinitializer.h"

/**
 * @brief A namespace to get a LevelsetReinitializer type based on a specified
 * constexpr.
 */
namespace LevelsetReinitializerSetup {

/**
 * @brief Function returning the typedef of a LevelsetReinitializer based on a
 * constexpr template.
 * @tparam LevelsetReinitializers The constexpr template parameter to specify
 * the exact LevelsetReinitializer type.
 */
template <LevelsetReinitializers> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<LevelsetReinitializers::Weno> {
  typedef WenoIterativeLevelsetReinitializer type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<LevelsetReinitializers::Min> {
  typedef MinIterativeLevelsetReinitializer type;
};

} // namespace LevelsetReinitializerSetup

#endif // LEVELSET_REINITIALIZER_SETUP_H
