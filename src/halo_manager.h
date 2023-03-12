//===------------------------ halo_manager.h ------------------------------===//
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
#ifndef HALO_MANAGER_H
#define HALO_MANAGER_H

#include "block_definitions/field_interface_definitions.h"
#include "boundary_condition/external_halo_manager.h"
#include "communication/communication_manager.h"
#include "communication/internal_halo_manager.h"

/**
 * @brief The halo manager class provides the functionality to handle halo
 * updates from one node to another (internal) or to update the halos lying in
 * external boundaries.
 */
class HaloManager {

private:
  Tree &tree_;
  ExternalHaloManager const &external_halo_manager_;
  InternalHaloManager &internal_halo_manager_;
  CommunicationManager const &communication_manager_;
  unsigned int const maximum_level_;

public:
  HaloManager() = delete;
  explicit HaloManager(Tree &tree,
                       ExternalHaloManager const &external_halo_manager,
                       InternalHaloManager &internal_halo_manager,
                       CommunicationManager const &communication_manager,
                       unsigned int const maximum_level);
  ~HaloManager() = default;
  HaloManager(HaloManager const &) = delete;
  HaloManager &operator=(HaloManager const &) = delete;
  HaloManager(HaloManager &&) = delete;
  HaloManager &operator=(HaloManager &&) = delete;

  void MaterialHaloUpdate(std::vector<unsigned int> const &levels_ascending,
                          MaterialFieldType const field_type,
                          bool const cut_jumps = false) const;
  void MaterialHaloUpdateOnLevel(unsigned int const level,
                                 MaterialFieldType const field_type,
                                 bool const cut_jumps = false) const;
  void MaterialHaloUpdateOnLmax(MaterialFieldType const field_type,
                                bool const cut_jumps = true) const;
  void MaterialHaloUpdateOnLmaxMultis(MaterialFieldType const field_type) const;

  void MaterialInternalHaloUpdateOnLevel(unsigned int const level,
                                         MaterialFieldType const field_type,
                                         bool const cut_jumps = false) const;
  void
  MaterialExternalHaloUpdateOnLevel(unsigned int const level,
                                    MaterialFieldType const field_type) const;

  /**
   * @brief Calls an interface tag halo update on Lmax only.
   * @tparam IDB Level-set buffer type.
   */
  template <InterfaceDescriptionBufferType IDB>
  void InterfaceTagHaloUpdateOnLmax() const {
    InterfaceTagHaloUpdateOnLevelList<IDB>({maximum_level_});
  }

  /**
   * @brief Adjusts the values in the interface tag buffer according to their
   * type. (symmetry, internal ...).
   * @param updated_levels The levels on which halos of nodes will be modified.
   * @tparam IDB Level-set buffer type.
   */
  template <InterfaceDescriptionBufferType IDB>
  void InterfaceTagHaloUpdateOnLevelList(
      std::vector<unsigned int> const &updated_levels) const {
    for (unsigned int const &level : updated_levels) {
      internal_halo_manager_.InterfaceTagHaloUpdateOnLevel(level, IDB);
      // Update of domain boundaries
      for (auto const &domain_boundary :
           communication_manager_.ExternalBoundaries(level)) {
        nid_t const id = std::get<0>(domain_boundary);
        BoundaryLocation const location = std::get<1>(domain_boundary);
        external_halo_manager_.UpdateInterfaceTagExternal(
            tree_.GetNodeWithId(id).GetInterfaceTags<IDB>(), location);
      }
    } // levels
  }

  void InterfaceHaloUpdateOnLevelList(
      std::vector<unsigned int> const updated_levels,
      InterfaceBlockBufferType const halo_type) const;
  void InterfaceHaloUpdateOnLmax(InterfaceBlockBufferType const type) const;
};

#endif // HALO_MANAGER_H
