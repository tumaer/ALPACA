//===-------------------- levelset_reinitializer.h ------------------------===//
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
#ifndef LEVELSET_REINITIALIZER_H
#define LEVELSET_REINITIALIZER_H

#include "halo_manager.h"
#include "levelset/geometry/geometry_calculator.h"
#include "levelset/geometry/geometry_calculator_marching_cubes.h"
#include "materials/material_manager.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief The class LevelsetReinitializer ensures the (signed-)distance property
 * of a level-set field.
 * @tparam DerivedLevelsetReinitializer Typename as template parameter due to
 * CRTP.
 */
template <typename DerivedLevelsetReinitializer> class LevelsetReinitializer {

  friend DerivedLevelsetReinitializer;

protected:
  HaloManager
      &halo_manager_; // TODO-19 NH Think about making it const (rats tail)
  LogWriter &logger_;

  /**
   * @brief The default constructor for a LevelsetReinitializer object.
   * @param halo_manager Instance to a HaloManager which provides MPI-related
   * methods.
   */
  explicit LevelsetReinitializer(HaloManager &halo_manager)
      : halo_manager_(halo_manager), logger_(LogWriter::Instance()) {
    // Empty Constructor, besides initializer list.
  }

public:
  LevelsetReinitializer() = delete;
  ~LevelsetReinitializer() = default;
  LevelsetReinitializer(LevelsetReinitializer const &) = delete;
  LevelsetReinitializer &operator=(LevelsetReinitializer const &) = delete;
  LevelsetReinitializer(LevelsetReinitializer &&) = delete;
  LevelsetReinitializer &operator=(LevelsetReinitializer &&) = delete;

  /**
   * @brief The method to reinitialize a level-set field.
   * @param nodes The nodes for which the level-set field is reinitialized.
   * @param levelset_type  Level-set field type that is reinitialized.
   * @param is_last_stage Flag if the method is called for the last stage of the
   * applied Runge-Kutta method.
   */
  void Reinitialize(std::vector<std::reference_wrapper<Node>> const &nodes,
                    InterfaceDescriptionBufferType const levelset_type,
                    bool const is_last_stage) const {
    static_cast<DerivedLevelsetReinitializer const &>(*this)
        .ReinitializeImplementation(nodes, levelset_type, is_last_stage);
  }
};

#endif // LEVELSET_REINITIALIZER_H
