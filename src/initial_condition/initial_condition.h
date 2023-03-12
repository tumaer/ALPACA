//===----------------------- initial_condition.h --------------------------===//
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
#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

#include <memory>

#include "levelset_initializer.h"
#include "materials/material_definitions.h"
#include "prime_state_initializer.h"
#include "topology/node_id_type.h"

/**
 * @brief The InitialCondition class is used to set the state of all cells
according to the user input at the beginning of the simulation. It serves as a
 *        simple handler deligating the work to the appropriate members.
´ */
class InitialCondition {

  std::unique_ptr<PrimeStateInitializer const> prime_state_initializer_;
  // for an arbitrary number of levelsets, this needs to be replaced by a vector
  // of unique pointers. Cannot be const due to non-const member functions.
  std::unique_ptr<LevelsetInitializer> levelset_initializer_;

public:
  InitialCondition() = delete;
  explicit InitialCondition(
      std::unique_ptr<PrimeStateInitializer const> prime_state_initializer,
      std::unique_ptr<LevelsetInitializer> levelset_initializer)
      : prime_state_initializer_(std::move(prime_state_initializer)),
        levelset_initializer_(std::move(levelset_initializer)) {
    /** Empty besides initializer list */
  }
  ~InitialCondition() = default;
  InitialCondition(InitialCondition const &) = delete;
  InitialCondition &operator=(InitialCondition const &) = delete;
  InitialCondition(InitialCondition &&) = delete;
  InitialCondition &operator=(InitialCondition &&) = delete;

  // Proxy functions to be called from outside
  void GetInitialPrimeStates(
      nid_t const node_id, MaterialName const material,
      double (&prime_state_buffer)[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()])
      const {
    prime_state_initializer_->GetInitialPrimeStates(node_id, material,
                                                    prime_state_buffer);
  };
  void GetInitialLevelset(
      nid_t const node_id,
      double (&levelset_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]) {
    levelset_initializer_->GetInitialLevelset(node_id, levelset_buffer);
  }
  std::vector<MaterialName> GetInitialMaterials(nid_t const node_id) {
    return levelset_initializer_->GetInitialMaterials(node_id);
  }
};

#endif // INITIAL_CONDITION_H
