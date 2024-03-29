//===----------------------- two_phase_manager.cpp ------------------------===//
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
#include "two_phase_manager.h"

#include "material_sign_capsule.h"
#include "utilities/buffer_operations_interface.h"

/**
 * @brief Default constructor for the TwoPhaseManager. Calls the default
 * constructor of the base class.
 * @param material_manager Instance of a material manager, which already has
 * been initialized according to the user input.
 * @param halo_manager Instance to a HaloManager which provides MPI-related
 * methods.
 */
TwoPhaseManager::TwoPhaseManager(MaterialManager const &material_manager,
                                 HaloManager &halo_manager)
    : MultiPhaseManager(material_manager, halo_manager) {
#ifndef PERFORMANCE
  // NH I think != 2 would block degeneration to single phase, but I'm not sure.
  if (material_manager.GetNumberOfMaterials() > 2) {
    throw std::logic_error(
        "Do not use TwoPhaseManager for more than two materials");
  }
#endif
}

/**
 * @brief Allow cut-cell mixing for a single level-set field satisfying the
 * signed-distance property. See also base class.
 * @param nodes See base class.
 */
void TwoPhaseManager::MixImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes) const {

  // TODO-19 JW: If integration is done on total cells, this halo update can
  // possibly be left out
  halo_manager_.MaterialHaloUpdateOnLmax(MaterialFieldType::Conservatives);

  for (Node &node : nodes) {
    cut_cell_mixer_.Mix(node);
    buffer_handler_.TransformToVolumeAveragedConservatives(node);
  }

  halo_manager_.MaterialHaloUpdateOnLmax(MaterialFieldType::Conservatives);

  /** Extension is done on the prime state buffer. Thus, cells in which we do
   * not extend have to be filled with prime states obtained from the integrated
   * conservatives from the right-hand side buffer. This leads to the following
   * occupation of the prime-state buffers: Real-material cells and cut-cells in
   * which we do not extend: Prime states obtained from integrated
   * conservatives. Extension-band cells and cut-cells in which we extend: Prime
   * states of the last RK stage.
   */
  for (Node &node : nodes) {
    buffer_handler_.CalculatePrimesFromIntegratedConservatives(node);
  }
}

/**
 * @brief Ensures a well-resolved single-level set field satisfying the
 * signed-distance property. Therefore, scale-separation and reinitialization
 * are performed.
 * @param nodes See base class.
 * @param is_last_stage See base class.
 */
void TwoPhaseManager::EnforceWellResolvedDistanceFunctionImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes,
    bool const is_last_stage) const {
  if (CC::ScaleSeparationActive() && is_last_stage) {
    scale_separator_.SeparateScales(
        nodes, InterfaceBlockBufferType::LevelsetReinitialized);
    /**
     * Since we also want to reinitialize scale-separated cells we cannot update
     * the interface tags at this place. We have to do it after the
     * reinitialization.
     */
  }

  bool const reinitialize =
      (ReinitializationConstants::ReinitializeOnlyInLastRkStage &&
       !is_last_stage)
          ? false
          : true;

  if (reinitialize) {
    levelset_reinitializer_.Reinitialize(
        nodes, InterfaceDescriptionBufferType::Reinitialized, is_last_stage);
  }

  if (is_last_stage) {
    for (Node &node : nodes) {
      buffer_handler_.AdaptConservativesToWellResolvedDistanceFunction(node);
    }
    halo_manager_.MaterialHaloUpdateOnLmax(MaterialFieldType::Conservatives);
    UpdateInterfaceTagsOnFinestLevel<
        InterfaceDescriptionBufferType::Reinitialized>(nodes);

    for (Node &node : nodes) {
      SetVolumeFractionBuffer<InterfaceDescriptionBufferType::Reinitialized>(
          node);
    }
    halo_manager_.InterfaceHaloUpdateOnLmax(
        InterfaceBlockBufferType::VolumeFractionReinitialized);
  }
}

/**
 * @brief Allows to extend material states to ghost cells.
 * @param nodes See base class.
 */
void TwoPhaseManager::ExtendPrimeStatesImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes) const {

  /** After the extension we have the following occupation of the prime-state
   * buffer: Real-material cells and cut-cells in which we do not extend: Prime
   * states obtained from integrated conservatives. Extension-band cells and
   * cut-cells in which we extend: Extended values. Other cells (Non-real
   * material or narrow band cells): 0.0. */
  ghost_fluid_extender_.Extend(nodes);

  /** After the extension the right-hand side buffer has to be populated with
   * the extended values. This leads to the following occupation of the
   * right-hand side buffer:
   * Real-fluid cells and cut-cells in which we do not extend: Integrated
   * conservatives. Extension-band cells and cut-cells in which we extend:
   * Conservatives obtained from the extended prime states. Other cells
   * (Non-real material or narrow band cells): 0.0. */
  for (Node &node : nodes) {
    buffer_handler_.CalculateConservativesFromExtendedPrimes(node);
  } // nodes

  /** To have integrated and extended conservatives also in the halo cells, we
   * have to do a Halo Update on the right-hand side buffer. */
  halo_manager_.MaterialHaloUpdateOnLmax(MaterialFieldType::Conservatives);
}

/**
 * @brief Extends interface parameter into narrow band.
 * @param nodes See base class.
 */
void TwoPhaseManager::ExtendInterfaceStatesImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes) const {

  // Calls the extension method
  // halo update of the interface quantities to extend (extender acts on
  // interface states)
  for (InterfaceState const is : IF::ISTE()) {
    halo_manager_.InterfaceHaloUpdateOnLmax(
        MapInterfaceStateToInterfaceBlockBufferType(is));
  }
  // Extend all interface states
  interface_extender_.Extend(nodes);
}

/**
 * @brief Calculates the volume fractions for the positive material on the inner
 * cells of a given node and saves them in the volume-fraction buffer. The
 * volume fractions are calculated using the reinitialized levelset buffer.
 * @param node The node for which the volume fraction buffer is set.
 * @tparam T The level set description for which the volume fraction buffer is
 * set.
 */
template <InterfaceDescriptionBufferType T>
void TwoPhaseManager::SetVolumeFractionBuffer(Node &node) const {

  std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceTags<T>();

  double const(&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock()
          .GetInterfaceDescriptionBuffer<T>()[InterfaceDescription::Levelset];
  double(&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock().GetInterfaceDescriptionBuffer<T>()
          [InterfaceDescription::VolumeFraction];

  for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        if (std::abs(interface_tags[i][j][k]) <=
            ITTI(IT::NewCutCell)) { // volume fraction calculation only makes
                                    // sense for old or new cut cells, for
                                    // non-cut cells it is trivial.
          volume_fraction[i][j][k] =
              geometry_calculator_.ComputeVolumeFraction(levelset, i, j, k);
        } else if (interface_tags[i][j][k] > ITTI(IT::NewCutCell)) {
          volume_fraction[i][j][k] = 1.0;
        } else {
          volume_fraction[i][j][k] = 0.0;
        }
      } // k
    }   // j
  }     // i
}

/**
 * See base class.
 * @param nodes See base class.
 * @param is_last_stage See base class.
 */
void TwoPhaseManager::UpdateIntegratedBufferImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes,
    bool const is_last_stage) const {

  /***
   * At this point in the algorithm time integration of the level-set field is
   * already done. To fully propagate the level-set field, interface tags and
   * volume fractions have to be updated. The interface tag update and the
   * calculation of volume fractions are based on the level-set values in the
   * integrated level-set buffer. Thus, after the advection of the level-set
   * field, the advected level-set values, which are currently in the right-hand
   * side level-set buffer, have to be copied to the integrated level-set
   * buffer.
   */
  BO::Interface::CopyInterfaceDescriptionBufferForNodeList<
      InterfaceDescriptionBufferType::RightHandSide,
      InterfaceDescriptionBufferType::Integrated>(nodes);
  UpdateInterfaceTagsOnFinestLevel<InterfaceDescriptionBufferType::Integrated>(
      nodes);

  if (CC::ScaleSeparationActive() && is_last_stage &&
      !ReinitializationConstants::ReinitializeAfterMixing) {
    scale_separator_.SeparateScales(
        nodes, InterfaceBlockBufferType::LevelsetIntegrated);
    /**
     * Since we also want to reinitialize scale-separated cells we cannot update
     * the interface tags at this place. We have to do it after the
     * reinitialization.
     */
  }

  bool reinitialize =
      ((ReinitializationConstants::ReinitializeOnlyInLastRkStage &&
        !is_last_stage) ||
       ReinitializationConstants::ReinitializeAfterMixing)
          ? false
          : true;
  if constexpr (CC::SolidBoundaryActive()) {
    for (Node &node : nodes) {
      for (auto &phase : node.GetPhases()) {
        if (material_manager_.IsSolidBoundary(phase.first)) {
          reinitialize = false;
        }
      }
    }
  }
  if (reinitialize) {
    levelset_reinitializer_.Reinitialize(
        nodes, InterfaceDescriptionBufferType::Integrated, is_last_stage);
  }

  UpdateInterfaceTagsOnFinestLevel<InterfaceDescriptionBufferType::Integrated>(
      nodes);

  // Set the volume fraction buffer according to the propagated level-set field.
  for (Node &node : nodes) {
    SetVolumeFractionBuffer<InterfaceDescriptionBufferType::Integrated>(node);
  }
  // A halo update for the volume fractions is necessary to have also correct
  // volume fractions in the halo cells.
  halo_manager_.InterfaceHaloUpdateOnLmax(
      InterfaceBlockBufferType::VolumeFractionIntegrated);
}

/**
 * See base class.
 * @param nodes See base class.
 */
void TwoPhaseManager::PropagateLevelsetImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes) const {

  BO::Interface::CopyInterfaceDescriptionBufferForNodeList<
      InterfaceDescriptionBufferType::Integrated,
      InterfaceDescriptionBufferType::Reinitialized,
      InterfaceDescription::Levelset>(nodes);
  BO::Interface::CopyInterfaceDescriptionBufferForNodeList<
      InterfaceDescriptionBufferType::Integrated,
      InterfaceDescriptionBufferType::Reinitialized,
      InterfaceDescription::VolumeFraction>(nodes);

  for (Node &node : nodes) {
    BO::CopySingleBuffer(
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Integrated>(),
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>());
  }
}

/**
 * @brief See base class.
 * @param nodes See base class.
 */
void TwoPhaseManager::InitializeVolumeFractionBufferImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes) const {
  for (Node &node : nodes) {
    SetVolumeFractionBuffer<InterfaceDescriptionBufferType::Reinitialized>(
        node);
  }
  halo_manager_.InterfaceHaloUpdateOnLmax(
      InterfaceBlockBufferType::VolumeFractionReinitialized);
}

/**
 * @brief See base class.
 * @param nodes See base class.
 */
void TwoPhaseManager::ObtainInterfaceStatesImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes,
    bool const reset_interface_states) const {
  if (reset_interface_states) {
    for (Node &node : nodes) {
      BO::SetFieldBuffer(node.GetInterfaceBlock().GetInterfaceStateBuffer(),
                         0.0);
    }
  }
  for (auto const &node : nodes) {
    // solve interface Riemann problem to obtain interface velocity and
    // interface exchange terms
    interface_state_calculator_.ObtainInterfaceStates(node);
  }

  for (InterfaceState const is : IF::ASOS()) {
    halo_manager_.InterfaceHaloUpdateOnLmax(
        MapInterfaceStateToInterfaceBlockBufferType(is));
  }

  if constexpr (InterfaceStateTreatmentConstants::ExtendInterfaceQuantities) {
    // extends interface quantities
    ExtendInterfaceStates(nodes);
  }
}
