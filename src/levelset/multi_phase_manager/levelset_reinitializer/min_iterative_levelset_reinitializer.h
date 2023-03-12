//===------------ min_iterative_levelset_reinitializer.h ------------------===//
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
#ifndef MIN_ITERATIVE_LEVELSET_REINITIALIZER_H
#define MIN_ITERATIVE_LEVELSET_REINITIALIZER_H

#include "iterative_levelset_reinitializer_base.h"

/**
 * @brief Provides functionality to reinitialize a level-set field as described
 * in \cite Sussman1994.
 */
class MinIterativeLevelsetReinitializer
    : public IterativeLevelsetReinitializerBase<
          MinIterativeLevelsetReinitializer> {

  friend IterativeLevelsetReinitializerBase;

  /**
   * The indicator whether subcell fix is used or not.
   */
  static constexpr bool subcell_fix_active_ = false;

protected:
  double ReinitializeSingleNodeImplementation(
      Node &node, InterfaceDescriptionBufferType const levelset_type,
      bool const) const;

public:
  MinIterativeLevelsetReinitializer() = delete;
  explicit MinIterativeLevelsetReinitializer(HaloManager &halo_manager);
  ~MinIterativeLevelsetReinitializer() = default;
  MinIterativeLevelsetReinitializer(MinIterativeLevelsetReinitializer const &) =
      delete;
  MinIterativeLevelsetReinitializer &
  operator=(MinIterativeLevelsetReinitializer const &) = delete;
  MinIterativeLevelsetReinitializer(MinIterativeLevelsetReinitializer &&) =
      delete;
  MinIterativeLevelsetReinitializer &
  operator=(MinIterativeLevelsetReinitializer &&) = delete;
};

#endif // MIN_ITERATIVE_LEVELSET_REINITIALIZER_H
