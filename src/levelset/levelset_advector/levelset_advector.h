//===----------------------- levelset_advector.h --------------------------===//
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
#ifndef LEVELSET_ADVECTOR_H
#define LEVELSET_ADVECTOR_H

#include "levelset/geometry/geometry_calculator_marching_cubes.h"
#include "topology/node.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief Provides functionality to propagate a level-set field in time by using
 * a predefined stencil.
 * @tparam DerivedLevelsetAdvector Typename as template parameter due to CRTP.
 */
template <typename DerivedLevelsetAdvector> class LevelsetAdvector {

  friend DerivedLevelsetAdvector;

public:
  // Private constructor only
  ~LevelsetAdvector() = default;
  LevelsetAdvector() = default;
  LevelsetAdvector(LevelsetAdvector const &) = delete;
  LevelsetAdvector &operator=(LevelsetAdvector const &) = delete;
  LevelsetAdvector(LevelsetAdvector &&) = delete;
  LevelsetAdvector &operator=(LevelsetAdvector &&) = delete;

  /**
   * @brief Calculates the right-hand side to solve the level-set advection
   * equation.
   * @param node The node for which the level-set field is propagated in time.
   */
  void Advect(Node &node) const {
    static_cast<DerivedLevelsetAdvector const &>(*this).AdvectImplementation(
        node);
  }
};

#endif // LEVELSET_ADVECTOR_H
