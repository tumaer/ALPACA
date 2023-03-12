//===----------------------- two_phase_manager.h --------------------------===//
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
#ifndef TWO_PHASE_MANAGER_H
#define TWO_PHASE_MANAGER_H

#include "buffer_handler.h"
#include "interface_tags/interface_tag_functions.h"
#include "multi_phase_manager.h"

/**
 * @brief The TwoPhaseManager provides functionality to perform two-phase flow
 * simulation by a single-level set method.
 */
class TwoPhaseManager : public MultiPhaseManager<TwoPhaseManager> {

  friend MultiPhaseManager;

private:
  template <InterfaceDescriptionBufferType T>
  void SetVolumeFractionBuffer(Node &node) const;

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

  /**
   * @brief Sets the interface tags on the finest level. Implementation for the
   * reinitialized buffer.
   * @param nodes_containing_level_set The nodes on the finest level, which have
   * a level-set block.
   * @tparam IDB The interface buffer type for the interface tag update.
   */
  template <InterfaceDescriptionBufferType IDB>
  void UpdateInterfaceTagsOnFinestLevel(
      std::vector<std::reference_wrapper<Node>> const
          &nodes_containing_level_set) const {

    for (Node &node : nodes_containing_level_set) {
      InterfaceTagFunctions::SetInternalCutCellTagsFromLevelset(
          node.GetInterfaceBlock().GetInterfaceDescriptionBuffer<IDB>()
              [InterfaceDescription::Levelset],
          node.GetInterfaceTags<IDB>());
    }

    halo_manager_.InterfaceTagHaloUpdateOnLmax<IDB>();

    for (Node &node : nodes_containing_level_set) {
      InterfaceTagFunctions::SetTotalInterfaceTagsFromCutCells(
          node.GetInterfaceTags<IDB>());
    }
    halo_manager_.InterfaceTagHaloUpdateOnLmax<IDB>();
  }

public:
  TwoPhaseManager() = delete;
  explicit TwoPhaseManager(MaterialManager const &material_manager,
                           HaloManager &halo_manager);
  ~TwoPhaseManager() = default;
  TwoPhaseManager(TwoPhaseManager const &) = delete;
  TwoPhaseManager &operator=(TwoPhaseManager const &) = delete;
  TwoPhaseManager(TwoPhaseManager &&) = delete;
  TwoPhaseManager &operator=(TwoPhaseManager &&) = delete;
};

#endif // TWO_PHASE_MANAGER_H
