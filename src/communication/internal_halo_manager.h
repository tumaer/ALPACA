//===--------------------- internal_halo_manager.h ------------------------===//
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
#ifndef INTERNAL_BOUNDARY_MANAGER_H
#define INTERNAL_BOUNDARY_MANAGER_H

#include "communication/communication_manager.h"
#include "communication/communication_statistics.h"

/**
 * @brief The InternalHaloManager is used within the domain, i.e. a classical
 * Halo. This class exchanges information of neighboring blocks by filling the
 * halos of the host with data from the domain of the neighbor. In the MR setup
 * two types of internal boundaries are possible, jump and no-jump; both are
 * treated by this class. In the first case communication with the parent of the
 * host node is needed in the latter communication with the direct neighbor is
 * needed. Works inter and intra MPI ranks.
 */
class InternalHaloManager {
private:
  Tree &tree_;
  TopologyManager &topology_;
  CommunicationManager
      &communication_manager_; // Cannot be const (for now NH TODO-19) because
                               // of new tagging system.
  unsigned int const number_of_materials_;

  // Helper function for the local halo filling of special buffers
  template <class T>
  inline void
  UpdateNoJumpLocal(T (&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                    T const (&partner_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                    BoundaryLocation const loc) const;

  template <class T>
  inline void
  ExtendClosestInternalValue(T (&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                             BoundaryLocation const loc) const;

  unsigned int UpdateMaterialJumpMpiSend(nid_t id,
                                         std::vector<MPI_Request> &requests,
                                         nid_t const remote_child_id,
                                         void *send_buffer,
                                         BoundaryLocation const loc,
                                         MaterialFieldType const filed_type);
  void UpdateMaterialJumpMpiRecv(nid_t id, std::vector<MPI_Request> &requests,
                                 BoundaryLocation const loc,
                                 MaterialFieldType const field_type);
  void UpdateMaterialJumpNoMpi(nid_t id, BoundaryLocation const loc,
                               MaterialFieldType const field_type);
  void UpdateMaterialHaloCellsMpiSend(nid_t id,
                                      std::vector<MPI_Request> &requests,
                                      BoundaryLocation const loc,
                                      MaterialFieldType const field_type);
  void UpdateMaterialHaloCellsMpiRecv(nid_t id,
                                      std::vector<MPI_Request> &requests,
                                      BoundaryLocation const loc,
                                      MaterialFieldType const field_type);
  void UpdateMaterialHaloCellsNoMpi(nid_t id, BoundaryLocation const loc,
                                    MaterialFieldType const field_type);

  void
  UpdateInterfaceHaloCellsMpiSend(nid_t id, std::vector<MPI_Request> &requests,
                                  InterfaceBlockBufferType const buffer_type,
                                  BoundaryLocation const loc);
  void
  UpdateInterfaceHaloCellsMpiRecv(nid_t id, std::vector<MPI_Request> &requests,
                                  InterfaceBlockBufferType const buffer_type,
                                  BoundaryLocation const loc);
  void UpdateInterfaceHaloCellsNoMpi(nid_t id,
                                     InterfaceBlockBufferType const buffer_type,
                                     BoundaryLocation const loc);

  void UpdateInterfaceTagHaloCellsMpiSend(
      nid_t id, std::vector<MPI_Request> &requests,
      InterfaceDescriptionBufferType const type, BoundaryLocation const loc);
  void UpdateInterfaceTagHaloCellsMpiRecv(
      nid_t id, std::vector<MPI_Request> &requests,
      InterfaceDescriptionBufferType const type, BoundaryLocation const loc);
  void UpdateInterfaceTagHaloCellsNoMpi(
      nid_t id, InterfaceDescriptionBufferType const buffer_type,
      BoundaryLocation const loc);

  void MpiMaterialHaloUpdateNoJump(
      std::vector<MPI_Request> &requests,
      std::vector<std::tuple<nid_t, BoundaryLocation,
                             InternalBoundaryType>> const &no_jump_boundaries,
      MaterialFieldType const field_type);
  void NoMpiMaterialHaloUpdate(
      std::vector<std::tuple<nid_t, BoundaryLocation,
                             InternalBoundaryType>> const &boundaries,
      MaterialFieldType const field_type);
  void MpiMaterialHaloUpdateJump(
      std::vector<MPI_Request> &requests,
      std::vector<std::tuple<nid_t, BoundaryLocation,
                             InternalBoundaryType>> const &no_jump_boundaries,
      std::vector<ExchangePlane> &jump_buffer_plane,
      std::vector<ExchangeStick> &jump_buffer_stick,
      std::vector<ExchangeCube> &jump_buffer_cube,
      MaterialFieldType const field_type);

  void NoMpiInterfaceTagHaloUpdate(
      std::vector<std::tuple<nid_t, BoundaryLocation,
                             InternalBoundaryType>> const &boundaries,
      InterfaceDescriptionBufferType const type);

  void MpiInterfaceTagHaloUpdate(
      std::vector<std::tuple<nid_t, BoundaryLocation,
                             InternalBoundaryType>> const &boundaries,
      InterfaceDescriptionBufferType const type,
      std::vector<MPI_Request> &requests);

  void NoMpiInterfaceHaloUpdate(
      std::vector<std::tuple<nid_t, BoundaryLocation,
                             InternalBoundaryType>> const &boundaries,
      InterfaceBlockBufferType const buffer_type);
  void MpiInterfaceHaloUpdate(
      std::vector<std::tuple<nid_t, BoundaryLocation,
                             InternalBoundaryType>> const &boundaries,
      InterfaceBlockBufferType const buffer_type,
      std::vector<MPI_Request> &requests);

public:
  InternalHaloManager() = delete;
  explicit InternalHaloManager(Tree &tree, TopologyManager &topology,
                               CommunicationManager &communication_manager,
                               unsigned int const number_of_materials);
  ~InternalHaloManager() = default;
  InternalHaloManager(InternalHaloManager const &) = delete;
  InternalHaloManager &operator=(InternalHaloManager const &) = delete;
  InternalHaloManager(InternalHaloManager &&) = delete;
  InternalHaloManager &operator=(InternalHaloManager &&) = delete;

  void MaterialHaloUpdateOnLevel(unsigned int const level,
                                 MaterialFieldType const field_type,
                                 bool const cut_jumps);

  void MaterialHaloUpdateOnMultis(MaterialFieldType const field_type);

  void InterfaceTagHaloUpdateOnLevel(unsigned int const level,
                                     InterfaceDescriptionBufferType const type);

  void InterfaceHaloUpdateOnLevel(unsigned int const level,
                                  InterfaceBlockBufferType const type);
};

#endif // INTERNAL_BOUNDARY_MANAGER_H
