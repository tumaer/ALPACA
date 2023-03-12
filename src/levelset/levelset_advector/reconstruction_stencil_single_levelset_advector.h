//===------- reconstruction_stencil_single_levelset_advector.h ------------===//
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
#ifndef RECONSTRUCTION_STENCIL_SINGLE_LEVELSET_ADVECTOR_H
#define RECONSTRUCTION_STENCIL_SINGLE_LEVELSET_ADVECTOR_H

#include "levelset_advector.h"

/**
 * @brief The ReconstructionStencilSingleLevelsetAdvector propagates a
 * single-level set field in time.
 */
class ReconstructionStencilSingleLevelsetAdvector
    : public LevelsetAdvector<ReconstructionStencilSingleLevelsetAdvector> {

  friend LevelsetAdvector;

  void AdvectImplementation(Node &node) const;

public:
  explicit ReconstructionStencilSingleLevelsetAdvector() = default;
  ~ReconstructionStencilSingleLevelsetAdvector() = default;
  ReconstructionStencilSingleLevelsetAdvector(
      ReconstructionStencilSingleLevelsetAdvector const &) = delete;
  ReconstructionStencilSingleLevelsetAdvector &
  operator=(ReconstructionStencilSingleLevelsetAdvector const &) = delete;
  ReconstructionStencilSingleLevelsetAdvector(
      ReconstructionStencilSingleLevelsetAdvector &&) = delete;
  ReconstructionStencilSingleLevelsetAdvector &
  operator=(ReconstructionStencilSingleLevelsetAdvector &&) = delete;
};

#endif // RECONSTRUCTION_STENCIL_SINGLE_LEVELSET_ADVECTOR_H
