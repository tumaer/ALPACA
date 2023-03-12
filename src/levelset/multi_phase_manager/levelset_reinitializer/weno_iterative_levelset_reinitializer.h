//===----------- weno_iterative_levelset_reinitializer.h ------------------===//
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
#ifndef WENO_ITERATIVE_LEVELSET_REINITIALIZER_H
#define WENO_ITERATIVE_LEVELSET_REINITIALIZER_H

#include "iterative_levelset_reinitializer_base.h"

/**
 * @brief Provides functionality to reinitialize a level-set field with a
 * reconstruction stencil.
 */
class WenoIterativeLevelsetReinitializer
    : public IterativeLevelsetReinitializerBase<
          WenoIterativeLevelsetReinitializer> {

  friend IterativeLevelsetReinitializerBase;

protected:
  double ReinitializeSingleNodeImplementation(
      Node &node, InterfaceDescriptionBufferType const levelset_type,
      bool const is_last_stage) const;

public:
  WenoIterativeLevelsetReinitializer() = delete;
  explicit WenoIterativeLevelsetReinitializer(HaloManager &halo_manager);
  ~WenoIterativeLevelsetReinitializer() = default;
  WenoIterativeLevelsetReinitializer(
      WenoIterativeLevelsetReinitializer const &) = delete;
  WenoIterativeLevelsetReinitializer &
  operator=(WenoIterativeLevelsetReinitializer const &) = delete;
  WenoIterativeLevelsetReinitializer(WenoIterativeLevelsetReinitializer &&) =
      delete;
  WenoIterativeLevelsetReinitializer &
  operator=(WenoIterativeLevelsetReinitializer &&) = delete;
};

#endif // WENO_ITERATIVE_LEVELSET_REINITIALIZER_H
