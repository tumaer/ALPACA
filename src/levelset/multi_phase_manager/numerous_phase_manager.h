//===------------------- numerous_phase_manager.h -------------------------===//
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
#ifndef NUMEROUS_PHASE_MANAGER_H
#define NUMEROUS_PHASE_MANAGER_H

#include "multi_phase_manager.h"

/**
 * @brief The NumerouesPhaseManager enables to simulate multi-phase flows with a
 * multi-level set method.
 */
class NumerousPhaseManager : public MultiPhaseManager<NumerousPhaseManager> {

  friend MultiPhaseManager;

  void MixImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes) const;
  void EnforceWellResolvedDistanceFunctionImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes,
      bool const is_last_stage = false) const;
  void ExtendPrimeStatesImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes) const;
  void ExtendInterfaceStatesImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes) const;
  void UpdateIntegratedBufferImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes,
      bool const is_last_stage) const;
  void PropagateLevelsetImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes) const;
  void InitializeVolumeFractionBufferImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes) const;
  void ObtainInterfaceStatesImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes,
      bool const reset_interface_states = false) const;

public:
  NumerousPhaseManager() = delete;
  explicit NumerousPhaseManager(MaterialManager const &material_manager,
                                HaloManager &halo_manager);
  ~NumerousPhaseManager() = default;
  NumerousPhaseManager(NumerousPhaseManager const &) = delete;
  NumerousPhaseManager &operator=(NumerousPhaseManager const &) = delete;
  NumerousPhaseManager(NumerousPhaseManager &&) = delete;
  NumerousPhaseManager &operator=(NumerousPhaseManager &&) = delete;
};

#endif // NUMEROUS_PHASE_MANAGER_H
