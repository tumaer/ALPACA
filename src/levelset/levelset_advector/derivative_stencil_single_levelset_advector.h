//===--------- derivative_stencil_single_levelset_advector.h --------------===//
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
#ifndef DERIVATIVE_STENCIL_SINGLE_LEVELSET_ADVECTOR_H
#define DERIVATIVE_STENCIL_SINGLE_LEVELSET_ADVECTOR_H

#include "levelset_advector.h"

/**
 * @brief The DerivativeStencilSingleLevelsetAdvector propagates a single-level
 * set field in time.
 */
class DerivativeStencilSingleLevelsetAdvector
    : public LevelsetAdvector<DerivativeStencilSingleLevelsetAdvector> {

  friend LevelsetAdvector;

  void AdvectImplementation(Node &node) const;

public:
  explicit DerivativeStencilSingleLevelsetAdvector() = default;
  ~DerivativeStencilSingleLevelsetAdvector() = default;
  DerivativeStencilSingleLevelsetAdvector(
      DerivativeStencilSingleLevelsetAdvector const &) = delete;
  DerivativeStencilSingleLevelsetAdvector &
  operator=(DerivativeStencilSingleLevelsetAdvector const &) = delete;
  DerivativeStencilSingleLevelsetAdvector(
      DerivativeStencilSingleLevelsetAdvector &&) = delete;
  DerivativeStencilSingleLevelsetAdvector &
  operator=(DerivativeStencilSingleLevelsetAdvector &&) = delete;
};

#endif // DERIVATIVE_STENCIL_SINGLE_LEVELSET_ADVECTOR_H
