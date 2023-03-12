//===----------------------- halo_manager.cpp -----------------------------===//
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
#include "halo_manager.h"
#include "communication/exchange_types.h"
#include "topology/id_information.h"

/**
 * @brief Default constructor for Halo Manager instance.
 * @param tree Instance to provide local (current rank) node information.
 * @param external_halo_manager Instance to provide halo updates on nodes having
 * an external boundary.
 * @param internal_halo_manager Instance to provide halo updates on nodes having
 * internal boundaries.
 * @param communication_manager Instance to access communication information.
 * @param maximum_level The maximum level used in the simulation.
 */
HaloManager::HaloManager(Tree &tree,
                         ExternalHaloManager const &external_halo_manager,
                         InternalHaloManager &internal_halo_manager,
                         CommunicationManager const &communication_manager,
                         unsigned int const maximum_level)
    : tree_(tree), external_halo_manager_(external_halo_manager),
      internal_halo_manager_(internal_halo_manager),
      communication_manager_(communication_manager),
      maximum_level_(maximum_level) {
  // Empty besides initializer list.
}

/**
 * @brief Adjusts the values in the halo cells, according to their type.
 * (symmetry, internal, ...).
 * @param levels_ascending The levels on which halos of nodes will be modified
 * in ascending order.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 * @param cut_jumps Decider if jump halos should be updated on all specified
 * level. If true: jumps will not be updated on the coarsest level in
 * "upddate_levels".
 */
void HaloManager::MaterialHaloUpdate(
    std::vector<unsigned int> const &levels_ascending,
    MaterialFieldType const field_type, bool const cut_jumps) const {
  std::vector<unsigned int> no_jump_update_levels(levels_ascending);
  /* NH 2017-02-20: It may be that no-jump halos are not to be updated on the
   * coarsest level in the input list. Therefore this level is handled
   * separately. Afterwards a normal Halo update is performed on all remaining
   * levels in the input list.
   */
  if (cut_jumps) {
    unsigned int no_jump_extra_level = no_jump_update_levels.front();
    no_jump_update_levels.erase(no_jump_update_levels.begin());
    MaterialHaloUpdateOnLevel(no_jump_extra_level, field_type, true);
  }
  for (unsigned int const level : no_jump_update_levels) {
    MaterialHaloUpdateOnLevel(level, field_type, false);
  }
}

/**
 * @brief Adjusts the material values in all halo cells, according to their
 * type.
 * @param level The level on which halos of nodes will be modified.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 * @param cut_jumps Decider if jump halos should be updated on specified level.
 * If true: jumps will not be updated on the current level.
 */
void HaloManager::MaterialHaloUpdateOnLevel(unsigned int const level,
                                            MaterialFieldType const field_type,
                                            bool const cut_jumps) const {
  MaterialInternalHaloUpdateOnLevel(level, field_type, cut_jumps);
  MaterialExternalHaloUpdateOnLevel(level, field_type);
}

/**
 * @brief Adjusts the values in internal halo cells, according to their type.
 * @param level The level on which halos of nodes will be modified.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 * @param cut_jumps Decider if jump halos should be updated on specified level.
 * If true: jumps will not be updated on the current level.
 */
void HaloManager::MaterialInternalHaloUpdateOnLevel(
    unsigned int const level, MaterialFieldType const field_type,
    bool const cut_jumps) const {
  internal_halo_manager_.MaterialHaloUpdateOnLevel(level, field_type,
                                                   cut_jumps);
}

/**
 * @brief Adjusts the material values in external halo cells, according to their
 * type.
 * @param level The level on which halos of nodes will be modified.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void HaloManager::MaterialExternalHaloUpdateOnLevel(
    unsigned int const level, MaterialFieldType const field_type) const {
  for (std::tuple<nid_t, BoundaryLocation> const &boundary :
       communication_manager_.ExternalBoundaries(level)) {
    external_halo_manager_.UpdateMaterialExternal(
        tree_.GetNodeWithId(std::get<0>(boundary)), field_type,
        std::get<1>(boundary));
  }
}

/**
 * @brief Adjusts the material values in the halo cells on the finest level,
 * according to their type.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 * @param cut_jumps Decider if jump halos should be updated on specified level.
 * If true: jumps will not be updated on the current level.
 * @note The default value for cut_jumps is true.
 */
void HaloManager::MaterialHaloUpdateOnLmax(MaterialFieldType const field_type,
                                           bool const cut_jumps) const {
  MaterialHaloUpdateOnLevel(maximum_level_, field_type, cut_jumps);
}

/**
 * @brief Adjusts the material values in the halo cells on the finest level for
 * multi-materials only, according to their type.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void HaloManager::MaterialHaloUpdateOnLmaxMultis(
    MaterialFieldType const field_type) const {
  internal_halo_manager_.MaterialHaloUpdateOnMultis(field_type);
  for (std::tuple<nid_t, BoundaryLocation> const &boundary :
       communication_manager_.ExternalMultiBoundaries()) {
    external_halo_manager_.UpdateMaterialExternal(
        tree_.GetNodeWithId(std::get<0>(boundary)), field_type,
        std::get<1>(boundary));
  }
}

/**
 * @brief Calls a interface halo update of the specified interface buffer on
 * Lmax (only!).
 * @param type The identifier of the buffer that is to be updated.
 */
void HaloManager::InterfaceHaloUpdateOnLmax(
    InterfaceBlockBufferType const type) const {
  // perform interface halo update on Lmax
  InterfaceHaloUpdateOnLevelList({maximum_level_}, type);
}

/**
 * @brief Adjusts the values in the stated interface block buffer according to
 * their type. (symmetry, internal ...).
 * @param updated_levels The levels on which halos of nodes will be modified.
 * @param type The identifier of the buffer that is to be updated.
 */
void HaloManager::InterfaceHaloUpdateOnLevelList(
    std::vector<unsigned int> const updated_levels,
    InterfaceBlockBufferType const type) const {
  for (auto const &level : updated_levels) {
    internal_halo_manager_.InterfaceHaloUpdateOnLevel(level, type);
    // Update of domain boundaries
    for (auto const &domain_boundary :
         communication_manager_.ExternalBoundaries(level)) {
      nid_t const id = std::get<0>(domain_boundary);
      BoundaryLocation const location = std::get<1>(domain_boundary);
      external_halo_manager_.UpdateLevelsetExternal(tree_.GetNodeWithId(id),
                                                    type, location);
    }
  } // levels
}
