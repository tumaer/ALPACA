//===--------------------- numerous_phase_manager.cpp ---------------------===//
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
#include "numerous_phase_manager.h"

/**
 * @brief Default constructor for the NumerousPhaseManager. Calls the default
 * constructor of the base class.
 * @param material_manager Instance of a material manager, which already has
 * been initialized according to the user input.
 * @param communicator Instance to a CommunicationManager which provides
 * MPI-related methods.
 */
NumerousPhaseManager::NumerousPhaseManager(
    MaterialManager const &material_manager, HaloManager &halo_manager)
    : MultiPhaseManager(material_manager, halo_manager) {
  // Empty Constructor, besides call of base class constructor.
}

/**
 * @brief Implements a cut-cell mixing procedure for multi-level set
 * simulations. See also base class.
 * @param nodes See base class.
 */
void NumerousPhaseManager::MixImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes) const {
  for (Node &node : nodes) {
    cut_cell_mixer_.Mix(node);
  }
}

/**
 * @brief See base class.
 * @param nodes See base class.
 * @param is_last_stage See base class.
 */
void NumerousPhaseManager::EnforceWellResolvedDistanceFunctionImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes,
    bool const is_last_stage) const {
  levelset_reinitializer_.Reinitialize(
      nodes, InterfaceDescriptionBufferType::Reinitialized, is_last_stage);
}

/**
 * @brief See base class.
 * @param node See base class.
 */
void NumerousPhaseManager::ExtendPrimeStatesImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes) const {
  ghost_fluid_extender_.Extend(nodes);
}

/**
 * @brief See base class.
 * @param node See base class.
 */
void NumerousPhaseManager::ExtendInterfaceStatesImplementation(
    std::vector<std::reference_wrapper<Node>> const &nodes) const {
  interface_extender_.Extend(nodes);
}

/**
 * @brief See base class.
 * @param nodes See base class.
 * @param stage See base class.
 */
void NumerousPhaseManager::UpdateIntegratedBufferImplementation(
    std::vector<std::reference_wrapper<Node>> const &,
    bool const is_last_stage) const {
  (void)is_last_stage;
  throw std::logic_error(
      "Not yet implemented: "
      "NumerousPhaseManager::UpdateIntegratedBufferImplementation");
}

/**
 * @brief See base class.
 * @param nodes See base class.
 * @param stage See base class.
 */
void NumerousPhaseManager::PropagateLevelsetImplementation(
    std::vector<std::reference_wrapper<Node>> const &) const {
  throw std::logic_error(
      "Not yet implemented: NumerousPhaseManager::PropagateLevelset");
}

/**
 * @brief See base class.
 * @param nodes See base class.
 */
void NumerousPhaseManager::InitializeVolumeFractionBufferImplementation(
    std::vector<std::reference_wrapper<Node>> const &) const {
  throw std::logic_error(
      "Not yet implemented: "
      "NumerousPhaseManager::InitializeVolumeFractionBuffer");
}

/**
 * @brief See base class.
 * @param nodes See base class.
 */
void NumerousPhaseManager::ObtainInterfaceStatesImplementation(
    std::vector<std::reference_wrapper<Node>> const &, bool const) const {
  throw std::logic_error(
      "Not yet implemented: NumerousPhaseManager::ObtainInterfaceStates");
}
