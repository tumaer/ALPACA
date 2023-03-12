//===-------------------- internal_halo_manager.cpp -----------------------===//
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
#include "communication/internal_halo_manager.h"
#include "communication/communication_manager.h"
#include "multiresolution/multiresolution.h"
#include "topology/id_information.h"

/**
 * @brief Standard Constructor for Internal Boundaries.
 * @param tree Instance providing information of local node arrangement.
 * @param topology Instance providing information of global node arrangement.
 * @param communication_manager Instance to carry out actual communication and
 * holding the internal and external neighbor/location relations.
 * @param number_of_materials Maximum number of materials in the simulation.
 */
InternalHaloManager::InternalHaloManager(
    Tree &tree, TopologyManager &topology,
    CommunicationManager &communication_manager,
    unsigned int const number_of_materials)
    : tree_(tree), topology_(topology),
      communication_manager_(communication_manager),
      number_of_materials_(number_of_materials) {
  // Empty besides initializer list
}

/**
 * @brief Adjusts the values in internal halo cells, according to their type.
 * @param level The level on which halos of nodes will be modified.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 * @param cut_jumps Decider if jump halos should be updated on specified level.
 * If true: jumps will not be updated on the current level.
 */
void InternalHaloManager::MaterialHaloUpdateOnLevel(
    unsigned int const level, MaterialFieldType const field_type,
    bool const cut_jumps) {
  std::vector<MPI_Request> requests;
  communication_manager_.GenerateNeighborRelationForHaloUpdate(level);

  // Non-Jump halo update
  // it is necessary that first the non-jump boundaries are carried out to
  // ensure that all parent nodes contain the correct information in their halo
  // cells
  MpiMaterialHaloUpdateNoJump(
      requests, communication_manager_.InternalBoundariesMpi(level),
      field_type);
  NoMpiMaterialHaloUpdate(communication_manager_.InternalBoundaries(level),
                          field_type);
  // Jump halo updates
  if (!cut_jumps) {
    // All direction-types need the same buffer size, but have different
    // DataTypes for sending the data, Conservatives_Plane_EW is also
    // representative for Conservatives_Plane_NS and Conservatives_Plane_TB and
    // so on. Proper cast is handled within SendJumpToChild in
    // internal_boundary_condition as of C++17 vector is guaranteed to be
    // contiguous in memory
    std::vector<ExchangePlane> jump_buffer_plane(
        MF::ANOF(field_type) * number_of_materials_ *
        communication_manager_.JumpSendCount(level, ExchangeType::Plane));
    std::vector<ExchangeStick> jump_buffer_stick(
        MF::ANOF(field_type) * number_of_materials_ *
        communication_manager_.JumpSendCount(level, ExchangeType::Stick));
    std::vector<ExchangeCube> jump_buffer_cube(
        MF::ANOF(field_type) * number_of_materials_ *
        communication_manager_.JumpSendCount(level, ExchangeType::Cube));
    MpiMaterialHaloUpdateJump(
        requests, communication_manager_.InternalBoundariesJumpMpi(level),
        jump_buffer_plane, jump_buffer_stick, jump_buffer_cube, field_type);
    NoMpiMaterialHaloUpdate(
        communication_manager_.InternalBoundariesJump(level), field_type);
    MPI_Waitall(
        requests.size(), requests.data(),
        MPI_STATUSES_IGNORE); // buffer-vectors need to be alive till this point
  } else {
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  }
  requests.clear();
}

void InternalHaloManager::MaterialHaloUpdateOnMultis(
    MaterialFieldType const field_type) {
  std::vector<MPI_Request> requests;
  communication_manager_.GenerateNeighborRelationForHaloUpdate(
      topology_.GetMaximumLevel());

  // update inner boundaries
  MpiMaterialHaloUpdateNoJump(
      requests, communication_manager_.InternalMultiBoundariesMpi(),
      field_type);
  NoMpiMaterialHaloUpdate(communication_manager_.InternalMultiBoundaries(),
                          field_type);
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  requests.clear();
}

/**
 * @brief Updates the interface tags in internal halo cells on the given level.
 * @param level Level on which the update is done.
 * @param type Level-set buffer type on which the update is done.
 */
void InternalHaloManager::InterfaceTagHaloUpdateOnLevel(
    unsigned int const level, InterfaceDescriptionBufferType const type) {
  std::vector<MPI_Request> requests;
  communication_manager_.GenerateNeighborRelationForHaloUpdate(level);
  // Non-Jump halo update
  // it is necessary that first the non-jump boundaries are carried out to
  // ensure that all parent nodes contain the correct information in their halo
  // cells
  MpiInterfaceTagHaloUpdate(communication_manager_.InternalBoundariesMpi(level),
                            type, requests);
  NoMpiInterfaceTagHaloUpdate(communication_manager_.InternalBoundaries(level),
                              type);
  // Jump halo update
  // (Inteface tag jumps are always handled locally, but might be in the mpi
  // buffer for MaterialHaloUpdates.)
  NoMpiInterfaceTagHaloUpdate(
      communication_manager_.InternalBoundariesJump(level), type);
  // Inteface tag jumps are always handled locally, but might be in the mpi
  // buffer for MaterialHaloUpdates.
  NoMpiInterfaceTagHaloUpdate(
      communication_manager_.InternalBoundariesJumpMpi(level), type);
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  requests.clear();
}

/**
 * @brief Updates the internal halo cells of a interface block buffer on the
 * given level.
 * @param level Level on which the update is done.
 * @param type The identifier of the buffer that is to be updated.
 */
void InternalHaloManager::InterfaceHaloUpdateOnLevel(
    unsigned int const level, InterfaceBlockBufferType const type) {
  std::vector<MPI_Request> requests;
  communication_manager_.GenerateNeighborRelationForHaloUpdate(level);
  // Non-Jump halo update
  // it is necessary that first the non-jump boundaries are carried out to
  // ensure that all parent nodes contain the correct information in their halo
  // cells
  MpiInterfaceHaloUpdate(communication_manager_.InternalBoundariesMpi(level),
                         type, requests);
  NoMpiInterfaceHaloUpdate(communication_manager_.InternalBoundaries(level),
                           type);
  // Jump halo update
  // (levelset jumps are always handled locally, but might be in the mpi buffer
  // for MaterialHaloUpdates.)
  NoMpiInterfaceHaloUpdate(communication_manager_.InternalBoundariesJump(level),
                           type);
  NoMpiInterfaceHaloUpdate(
      communication_manager_.InternalBoundariesJumpMpi(level), type);
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  requests.clear();
}

/**
 * @brief Updates Jump Halo Cells, called on child, if parent is on another MPI
 * rank.
 * @param id The id  of the node to be updated.
 * @param requests Container to store the requests opened by this function.
 * Indirect return parameter.
 * @param loc BoundaryLocation to be updated.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void InternalHaloManager::UpdateMaterialJumpMpiRecv(
    nid_t const id, std::vector<MPI_Request> &requests,
    BoundaryLocation const loc, MaterialFieldType const field_type) {
  Node &node = tree_.GetNodeWithId(id);
  // parent always contains material of child
  for (auto const material : topology_.GetMaterialsOfNode(id)) {
    Block &host_block = node.GetPhaseByMaterial(material);

    nid_t const parent_id = ParentIdOfNode(id);
    int const sender_rank = topology_.GetRankOfNode(parent_id);

    // MPI
    MPI_Datatype recv_type =
        communication_manager_.RecvDatatype(loc, DatatypeForMpi::Double);
    switch (field_type) {
    case MaterialFieldType::Conservatives: {
      communication_manager_.Recv(&host_block.GetRightHandSideBuffer(),
                                  MF::ANOE(), recv_type, sender_rank, requests);
    } break;
    case MaterialFieldType::PrimeStates: {
      communication_manager_.Recv(&host_block.GetPrimeStateBuffer(), MF::ANOP(),
                                  recv_type, sender_rank, requests);
    } break;
#ifndef PERFORMANCE
    case MaterialFieldType::Parameters: {
      communication_manager_.Recv(&host_block.GetParameterBuffer(), MF::ANOPA(),
                                  recv_type, sender_rank, requests);
    } break;
    default:
      throw std::logic_error("Material field type not known!");
#else
    default: /* MaterialFieldType::Parameters: */
      communication_manager_.Recv(&host_block.GetParameterBuffer(), MF::ANOPA(),
                                  recv_type, sender_rank, requests);
#endif
    }
  }
}

/**
 * @brief Updates Jump Halo Cells, called on child if parent and child are on
 * the same MPI rank.
 * @param id The id of the node to be updated.
 * @param loc BoundaryLocation to be updated.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void InternalHaloManager::UpdateMaterialJumpNoMpi(
    nid_t const id, BoundaryLocation const loc,
    MaterialFieldType const field_type) {
  Node &node = tree_.GetNodeWithId(id);
  // parent always contains material of child
  for (auto const material : topology_.GetMaterialsOfNode(id)) {
    Block &host_block = node.GetPhaseByMaterial(material);

    nid_t const parent_id = ParentIdOfNode(id);
    auto const start_indices =
        communication_manager_.GetStartIndicesHaloRecv(loc);
    auto const halo_size = communication_manager_.GetHaloSize(loc);

    Block const &parent_block =
        tree_.GetNodeWithId(parent_id).GetPhaseByMaterial(material);

    switch (field_type) {
    case MaterialFieldType::Conservatives: {
      for (Equation const e : MF::ASOE()) {
        Multiresolution::Prediction(parent_block.GetRightHandSideBuffer(e),
                                    host_block.GetRightHandSideBuffer(e), id,
                                    start_indices[0], halo_size[0],
                                    start_indices[1], halo_size[1],
                                    start_indices[2], halo_size[2]);
      }
    } break;
    case MaterialFieldType::PrimeStates: {
      for (PrimeState const ps : MF::ASOP()) {
        Multiresolution::Prediction(parent_block.GetPrimeStateBuffer(ps),
                                    host_block.GetPrimeStateBuffer(ps), id,
                                    start_indices[0], halo_size[0],
                                    start_indices[1], halo_size[1],
                                    start_indices[2], halo_size[2]);
      }
    } break;
#ifndef PERFORMANCE
    case MaterialFieldType::Parameters: {
      for (Parameter const pa : MF::ASOPA()) {
        Multiresolution::Prediction(parent_block.GetParameterBuffer(pa),
                                    host_block.GetParameterBuffer(pa), id,
                                    start_indices[0], halo_size[0],
                                    start_indices[1], halo_size[1],
                                    start_indices[2], halo_size[2]);
      }
    } break;
    default:
      throw std::logic_error("Material field type not known!");
#else
    default: /* MaterialFieldType::Parameters: */ {
      for (Parameter const pa : MF::ASOPA()) {
        Multiresolution::Prediction(parent_block.GetParameterBuffer(pa),
                                    host_block.GetParameterBuffer(pa), id,
                                    start_indices[0], halo_size[0],
                                    start_indices[1], halo_size[1],
                                    start_indices[2], halo_size[2]);
      }
    }
#endif
    }
  }
}

/**
 * @brief Set the state of the Halo cells according to the Boundary Type.
 * @param id The id of the node to be updated.
 * @param loc BoundaryLocation to be updated.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void InternalHaloManager::UpdateMaterialHaloCellsNoMpi(
    nid_t const id, BoundaryLocation const loc,
    MaterialFieldType const field_type) {
  Node &node = tree_.GetNodeWithId(id);
  unsigned int const number_of_fields = MF::ANOF(field_type);
  for (auto const material : topology_.GetMaterialsOfNode(id)) {
    nid_t const neighbor_id = topology_.GetTopologyNeighborId(id, loc);
    if (topology_.NodeContainsMaterial(neighbor_id, material)) {
      Block &host_block = node.GetPhaseByMaterial(material);

      Block const &partner_block =
          tree_.GetNodeWithId(neighbor_id).GetPhaseByMaterial(material);

      for (unsigned int field_index = 0; field_index < number_of_fields;
           ++field_index) {
        double(&host_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] =
            host_block.GetFieldBuffer(field_type, field_index);
        double const(&partner_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] =
            partner_block.GetFieldBuffer(field_type, field_index);
        UpdateNoJumpLocal(host_cells, partner_cells, loc);
      }
    }
  }
}

/**
 * @brief Set the state of the Halo cells according to the Boundary Type. Only
 * Updates via MPI are executed by this function.
 * @param id The id of the node to be updated.
 * @param requests In case the Boundary Type requires MPI Communication
 * asynchronous communication is supported via a reference to a MPI_Request
 * vector.
 * @param loc BoundaryLocation to be updated.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void InternalHaloManager::UpdateMaterialHaloCellsMpiSend(
    nid_t const id, std::vector<MPI_Request> &requests,
    BoundaryLocation const loc, MaterialFieldType const field_type) {
  Node &node = tree_.GetNodeWithId(id);
  nid_t const neighbor_id = topology_.GetTopologyNeighborId(id, loc);
  int const rank_of_neighbor = topology_.GetRankOfNode(neighbor_id);

  for (auto const material : topology_.GetMaterialsOfNode(id)) {
    if (topology_.NodeContainsMaterial(neighbor_id, material)) {
      Block const &host_block = node.GetPhaseByMaterial(material);
      switch (field_type) {
      case MaterialFieldType::Conservatives: {
        MPI_Datatype send_type =
            communication_manager_.SendDatatype(loc, DatatypeForMpi::Double);
        communication_manager_.Send(&host_block.GetRightHandSideBuffer(),
                                    MF::ANOE(), send_type, rank_of_neighbor,
                                    requests);
      } break;
      case MaterialFieldType::PrimeStates: {
        MPI_Datatype send_type =
            communication_manager_.SendDatatype(loc, DatatypeForMpi::Double);
        communication_manager_.Send(&host_block.GetPrimeStateBuffer(),
                                    MF::ANOP(), send_type, rank_of_neighbor,
                                    requests);
      } break;
#ifndef PERFORMANCE
      case MaterialFieldType::Parameters: {
        MPI_Datatype send_type =
            communication_manager_.SendDatatype(loc, DatatypeForMpi::Double);
        communication_manager_.Send(&host_block.GetParameterBuffer(),
                                    MF::ANOPA(), send_type, rank_of_neighbor,
                                    requests);
      } break;
      default:
        throw std::logic_error("Material field type not known!");
#else
      default: /* MaterialFieldType::Parameters: */ {
        MPI_Datatype send_type =
            communication_manager_.SendDatatype(loc, DatatypeForMpi::Double);
        communication_manager_.Send(&host_block.GetParameterBuffer(),
                                    MF::ANOPA(), send_type, rank_of_neighbor,
                                    requests);
      }
#endif
      }
    }
  }
}

/**
 * @brief Set the state of the Halo cells according to the Boundary Type. Only
 * Updates via MPI are executed by this function.
 * @param id The id of the node to be updated.
 * @param requests In case the Boundary Type requires MPI Communication
 * asynchronous communication is supported via a reference to a MPI_Request
 * vector.
 * @param loc BoundaryLocation to be updated.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void InternalHaloManager::UpdateMaterialHaloCellsMpiRecv(
    nid_t const id, std::vector<MPI_Request> &requests,
    BoundaryLocation const loc, MaterialFieldType const field_type) {
  Node &node = tree_.GetNodeWithId(id);
  nid_t const neighbor_id = topology_.GetTopologyNeighborId(id, loc);
  int const rank_of_neighbor = topology_.GetRankOfNode(neighbor_id);

  for (auto const material : topology_.GetMaterialsOfNode(id)) {
    if (topology_.NodeContainsMaterial(neighbor_id, material)) {
      Block &host_block = node.GetPhaseByMaterial(material);
      MPI_Datatype recv_type =
          communication_manager_.RecvDatatype(loc, DatatypeForMpi::Double);
      switch (field_type) {
      case MaterialFieldType::Conservatives: {
        communication_manager_.Recv(&host_block.GetRightHandSideBuffer(),
                                    MF::ANOE(), recv_type, rank_of_neighbor,
                                    requests);
      } break;
      case MaterialFieldType::PrimeStates: {
        communication_manager_.Recv(&host_block.GetPrimeStateBuffer(),
                                    MF::ANOP(), recv_type, rank_of_neighbor,
                                    requests);
      } break;
      case MaterialFieldType::Parameters: {
        communication_manager_.Recv(&host_block.GetParameterBuffer(),
                                    MF::ANOPA(), recv_type, rank_of_neighbor,
                                    requests);
      } break;
      default:
        throw std::logic_error("Material field type not known!");
      }
    }
  }
}

/**
 * @brief Updates a no-jump boundary with both nodes present on this rank (thus
 * local).
 * @param host_buffer Reference to the buffer that is to be updated.
 * @param partner_buffer Reference to the buffer the data should be taken from.
 * @tparam T Base type of the buffers.
 */
template <class T>
inline void InternalHaloManager::UpdateNoJumpLocal(
    T (&__restrict host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
    T const (&__restrict partner_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
    BoundaryLocation const loc) const {
  auto const send_indices =
      communication_manager_.GetStartIndicesHaloSend(OppositeDirection(loc));
  auto const recv_indices = communication_manager_.GetStartIndicesHaloRecv(loc);
  auto const size = communication_manager_.GetHaloSize(loc);

  // NH Int as rat's tail form MPI which does not allow uint.
  for (int i = 0; i < size[0]; ++i) {
    for (int j = 0; j < size[1]; ++j) {
      for (int k = 0; k < size[2]; ++k) {
        host_buffer[i + recv_indices[0]][j + recv_indices[1]]
                   [k + recv_indices[2]] =
                       partner_buffer[i + send_indices[0]][j + send_indices[1]]
                                     [k + send_indices[2]];
      }
    }
  }
}

/**
 * @brief Fills the halo cells at the specified location with internal value
 * closest to the respective cell.
 * @param loc BoundaryLocation to be updated.
 * @tparam T Base type of the buffers.
 */
template <class T>
inline void InternalHaloManager::ExtendClosestInternalValue(
    T (&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
    BoundaryLocation const loc) const {
  auto const start_indices =
      communication_manager_.GetStartIndicesHaloRecv(loc);
  auto const halo_size = communication_manager_.GetHaloSize(loc);

  // NH int as rat's tail form MPI which does not allow uint.
  for (unsigned int i = start_indices[0];
       i < static_cast<unsigned int>(start_indices[0] + halo_size[0]); ++i) {
    for (unsigned int j = start_indices[1];
         j < static_cast<unsigned int>(start_indices[1] + halo_size[1]); ++j) {
      for (unsigned int k = start_indices[2];
           k < static_cast<unsigned int>(start_indices[2] + halo_size[2]);
           ++k) {
        unsigned int ii =
            i < CC::FICX() ? CC::FICX() : (i > CC::LICX() ? CC::LICX() : i);
        unsigned int jj =
            j < CC::FICY() ? CC::FICY() : (j > CC::LICY() ? CC::LICY() : j);
        unsigned int kk =
            k < CC::FICZ() ? CC::FICZ() : (k > CC::LICZ() ? CC::LICZ() : k);
        host_buffer[i][j][k] = host_buffer[ii][jj][kk];
      }
    }
  }
}

/**
 * @brief Sends the state of the Halo cells according to the Boundary Type.
 * @param id The id of the node to be updated.
 * @param requests Asynchronous communication is supported via a reference to a
 * MPI_Request vector.
 * @param buffer_type Type of buffer on the levelset block that should be
 * updated.
 * @param loc BoundaryLocation to be updated.
 */
void InternalHaloManager::UpdateInterfaceHaloCellsMpiSend(
    nid_t const id, std::vector<MPI_Request> &requests,
    InterfaceBlockBufferType const buffer_type, BoundaryLocation const loc) {
  Node &node = tree_.GetNodeWithId(id);
  nid_t const host_id = id;
  nid_t const neighbor_id = topology_.GetTopologyNeighborId(host_id, loc);
  /*  NH TODO-19 "node.HasLevelset() if" should be avoided, therefore different
   * neighbor relations in CommunicationManger needed for levelset vs.
   * Material/Tag Halo updates.
   */
  if (topology_.IsNodeMultiPhase(neighbor_id) && node.HasLevelset()) {
    int const rank_of_neighbor = topology_.GetRankOfNode(neighbor_id);
    double const(&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetBuffer(buffer_type);
    MPI_Datatype send_type =
        communication_manager_.SendDatatype(loc, DatatypeForMpi::Double);
    communication_manager_.Send(host_buffer, 1, send_type, rank_of_neighbor,
                                requests);
  }
}

/**
 * @brief Receives the state of the interface halo cells according to the
 * Boundary Type.
 * @param id The id of the node to be updated.
 * @param requests Asynchronous communication is supported via a reference to a
 * MPI_Request vector.
 * @param buffer_type Type of buffer on the levelset block that should be
 * updated.
 * @param loc BoundaryLocation to be updated.
 */
void InternalHaloManager::UpdateInterfaceHaloCellsMpiRecv(
    nid_t const id, std::vector<MPI_Request> &requests,
    InterfaceBlockBufferType const buffer_type, BoundaryLocation const loc) {
  Node &node = tree_.GetNodeWithId(id);
  nid_t const host_id = id;
  nid_t const neighbor_id = topology_.GetTopologyNeighborId(host_id, loc);
  /*  NH TODO-19 "node.HasLevelset() if" should be avoided, therefore different
   * neighbor relations in CommunicationManger needed for levelset vs.
   * Material/Tag Halo updates.
   */
  if (node.HasLevelset()) {
    if (topology_.IsNodeMultiPhase(neighbor_id)) {
      int const rank_of_neighbor = topology_.GetRankOfNode(neighbor_id);
      double(&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetBuffer(buffer_type);
      MPI_Datatype recv_type =
          communication_manager_.RecvDatatype(loc, DatatypeForMpi::Double);
      communication_manager_.Recv(host_buffer, 1, recv_type, rank_of_neighbor,
                                  requests);
    } else {
      double(&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetBuffer(buffer_type);
      ExtendClosestInternalValue(host_buffer, loc);
    }
  }
}

/**
 * @brief Performs all Halo Updates that can be done without communication.
 * @param id The id of the node to be updated.
 * @param buffer_type Type of buffer on the levelset block that should be
 * updated.
 * @param loc BoundaryLocation to be updated.
 */
void InternalHaloManager::UpdateInterfaceHaloCellsNoMpi(
    nid_t const id, InterfaceBlockBufferType const buffer_type,
    BoundaryLocation const loc) {
  Node &node = tree_.GetNodeWithId(id);
  nid_t const host_id = id;
  nid_t const neighbor_id = topology_.GetTopologyNeighborId(host_id, loc);

  /*  NH TODO-19 "node.HasLevelset() if" should be avoided, therefore different
   * neighbor relations in CommunicationManger needed for levelset vs.
   * Material/Tag Halo updates.
   */
  if (node.HasLevelset()) {
    if (topology_.NodeExists(neighbor_id) &&
        topology_.IsNodeMultiPhase(neighbor_id)) {
      double(&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetBuffer(buffer_type);
      Node const &neighbor = tree_.GetNodeWithId(neighbor_id);
      if (neighbor.HasLevelset()) {
        double const(&partner_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
            neighbor.GetInterfaceBlock().GetBuffer(buffer_type);
        UpdateNoJumpLocal(host_buffer, partner_buffer, loc);
      } else {
        ExtendClosestInternalValue(host_buffer, loc);
      }
    } else {
      double(&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetBuffer(buffer_type);
      ExtendClosestInternalValue(host_buffer, loc);
    }
  }
}

/**
 * @brief Sends the state of the interface tags according to the Boundary Type.
 * @param id_node The node and its id to be updated.
 * @param requests Asynchronous communication is supported via a reference to a
 * MPI_Request vector.
 * @param type Level-set buffer type on which the update is done.
 * @param loc BoundaryLocation to be updated.
 */
void InternalHaloManager::UpdateInterfaceTagHaloCellsMpiSend(
    nid_t const id, std::vector<MPI_Request> &requests,
    InterfaceDescriptionBufferType const type, BoundaryLocation const loc) {
  Node &node = tree_.GetNodeWithId(id);
  nid_t const host_id = id;
  nid_t const neighbor_id = topology_.GetTopologyNeighborId(host_id, loc);

  int const rank_of_neighbor = topology_.GetRankOfNode(neighbor_id);
  std::int8_t const(&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceTags(type);
  MPI_Datatype send_type =
      communication_manager_.SendDatatype(loc, DatatypeForMpi::Byte);
  communication_manager_.Send(host_buffer, 1, send_type, rank_of_neighbor,
                              requests);
}

/**
 * @brief Receives the state of the interface halo cells according to the
 * Boundary Type.
 * @param id The id of the node to be updated.
 * @param requests Asynchronous communication is supported via a reference to a
 * MPI_Request vector.
 * @param type Level-set buffer type on which the update is done.
 * @param loc BoundaryLocation to be updated.
 */
void InternalHaloManager::UpdateInterfaceTagHaloCellsMpiRecv(
    nid_t const id, std::vector<MPI_Request> &requests,
    InterfaceDescriptionBufferType const type, BoundaryLocation const loc) {
  Node &node = tree_.GetNodeWithId(id);
  nid_t const host_id = id;
  nid_t const neighbor_id = topology_.GetTopologyNeighborId(host_id, loc);

  int const rank_of_neighbor = topology_.GetRankOfNode(neighbor_id);

  std::int8_t(&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceTags(type);
  MPI_Datatype recv_type =
      communication_manager_.RecvDatatype(loc, DatatypeForMpi::Byte);
  communication_manager_.Recv(host_buffer, 1, recv_type, rank_of_neighbor,
                              requests);
}

/**
 * @brief Performs all interface halo updates that can be done without
 * communication.
 * @param id The id of the node to be updated.
 * @param type Level-set buffer type on which the update is done.
 * @param loc BoundaryLocation to be updated.
 */
void InternalHaloManager::UpdateInterfaceTagHaloCellsNoMpi(
    nid_t const id, InterfaceDescriptionBufferType const type,
    BoundaryLocation const loc) {
  Node &node = tree_.GetNodeWithId(id);
  nid_t const host_id = id;
  nid_t const neighbor_id = topology_.GetTopologyNeighborId(host_id, loc);

  if (topology_.NodeExists(neighbor_id)) {
    std::int8_t(&host_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags(type);
    std::int8_t const(&partner_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        tree_.GetNodeWithId(neighbor_id).GetInterfaceTags(type);
    UpdateNoJumpLocal(host_buffer, partner_buffer, loc);
  } else { // In case the neighbour does not exist, e.g. at jumps, the closest
           // internal value (bulk phase) has to be extended.
    ExtendClosestInternalValue(node.GetInterfaceTags(type), loc);
  }
}

/**
 * @brief Communication partner for UpdateMaterialJumpMpiRecv in case of a jump
 * boundary.
 * @param id The id of the node to be updated.
 * @param requests Asynchronous communication is supported via a reference to a
 * MPI_Request vector.
 * @param remote_child_id Id of the remote child.
 * @param send_buffer necessary for asynchronous communication, as buffer to
 * store temporary data. It has the size of Conservatives_Slice_EW but might
 * also be used for Conservatives_Slice_NS, Conservatives_Slice_TB.
 * @param loc BoundaryLocation to be updated.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
unsigned int InternalHaloManager::UpdateMaterialJumpMpiSend(
    nid_t const id, std::vector<MPI_Request> &requests,
    nid_t const remote_child_id, void *send_buffer, BoundaryLocation const loc,
    MaterialFieldType const field_type) {
  Node &node = tree_.GetNodeWithId(id);

  // NH ints needed here as a rat's tail of MPI, where uint is not allowed.
  auto const start_indices =
      communication_manager_.GetStartIndicesHaloRecv(loc);
  auto const halo_size = communication_manager_.GetHaloSize(loc);
  unsigned int const number_of_fields = MF::ANOF(field_type);
  int const single_block_size = halo_size[0] * halo_size[1] * halo_size[2];
  int const whole_block_size = number_of_fields * single_block_size;

  int material_number = 0;
  int const child_rank = topology_.GetRankOfNode(remote_child_id);

  for (MaterialName const material : topology_.GetMaterialsOfNode(id)) {

    if (topology_.NodeContainsMaterial(remote_child_id, material)) {
      Block const &parent_block = node.GetPhaseByMaterial(material);

      double *send_buffer_double = (double *)send_buffer;

      MPI_Datatype const send_type =
          communication_manager_.JumpPlaneSendDatatype(loc);
      int block_pos = material_number * whole_block_size;

      double child_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
      for (unsigned int field_index = 0; field_index < number_of_fields;
           ++field_index) {
        Multiresolution::Prediction(
            parent_block.GetFieldBuffer(field_type, field_index), child_buffer,
            remote_child_id, start_indices[0], halo_size[0], start_indices[1],
            halo_size[1], start_indices[2], halo_size[2]);
        /* in order to buffer only valid values, we need to copy the values into
         * the send buffer, this operation is dependent on the direction of the
         * jump the static_cast is used to specify the datatype of the memory
         * necessary to store the current temporary results. As we have to
         * distinguish the different directions, different types are used, but
         * all have the same total amount of bytes
         */
        for (int i = 0; i < halo_size[0]; i++) {
          for (int j = 0; j < halo_size[1]; j++) {
            for (int k = 0; k < halo_size[2]; k++) {
              int const cell_pos = material_number * whole_block_size +
                                   field_index * single_block_size +
                                   i * halo_size[1] * halo_size[2] +
                                   j * halo_size[2] + k;
              send_buffer_double[cell_pos] =
                  child_buffer[start_indices[0] + i][start_indices[1] + j]
                              [start_indices[2] + k];
            }
          }
        }
      }
      material_number++;
      communication_manager_.Send(&send_buffer_double[block_pos],
                                  number_of_fields, send_type, child_rank,
                                  requests);
    }
  }
  // return increment of field buffers sent
  return material_number * number_of_fields;
}

/**
 * @brief Method used to execute the Halo Update for all jump boundaries with
 * MPI Communication.
 * @param requests vector of communication request (handle), new handle will be
 * added at the end of the vector.
 * @param boundaries Description of internal boundaries with jump.
 * @param jump_buffer_plane Vector containing all buffers used for the jump
 * exchange of plane boundaries.
 * @param jump_buffer_stick Vector containing all buffers used for the jump
 * exchange of stick boundaries.
 * @param jump_buffer_cube Vector containing all buffers used for the jump
 * exchange of cube boundaries.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void InternalHaloManager::MpiMaterialHaloUpdateJump(
    std::vector<MPI_Request> &requests,
    std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const
        &boundaries,
    std::vector<ExchangePlane> &jump_buffer_plane,
    std::vector<ExchangeStick> &jump_buffer_stick,
    std::vector<ExchangeCube> &jump_buffer_cube,
    MaterialFieldType const field_type) {
  unsigned int jump_send_counter_plane = 0;
  unsigned int jump_send_counter_stick = 0;
  unsigned int jump_send_counter_cube = 0;

  for (auto const &boundary : boundaries) {
    nid_t const id = std::get<0>(boundary);
    BoundaryLocation const location = std::get<1>(boundary);

    switch (std::get<2>(boundary)) {
    case InternalBoundaryType::JumpBoundaryMpiRecv: {
      CommunicationStatistics::jump_halos_recv_++;
      UpdateMaterialJumpMpiRecv(id, requests, location, field_type);
    } break;
#ifndef PERFORMANCE
    case InternalBoundaryType::JumpBoundaryMpiSend: {
      CommunicationStatistics::jump_halos_send_++;
      nid_t const parent_id = ParentIdOfNode(id);

      // buffer size is dependent on jump location
      if (LTI(location) <= LTI(BoundaryLocation::Bottom)) {
        // Plane
        jump_send_counter_plane += UpdateMaterialJumpMpiSend(
            parent_id, requests, id,
            &jump_buffer_plane[jump_send_counter_plane], location, field_type);
      } else if (LTI(location) <= LTI(BoundaryLocation::SouthWest)) {
        // Stick
        jump_send_counter_stick += UpdateMaterialJumpMpiSend(
            parent_id, requests, id,
            &jump_buffer_stick[jump_send_counter_stick], location, field_type);
      } else {
        // Cube
        jump_send_counter_cube += UpdateMaterialJumpMpiSend(
            parent_id, requests, id, &jump_buffer_cube[jump_send_counter_cube],
            location, field_type);
      }
    } break;
    default:
      throw std::logic_error(" Jump Halo update: boundaryType does not need "
                             "MPI, is noJump or is unknown");
#else
    default: /* InternalBoundaryType::JumpBoundaryMpiSend */ {
      CommunicationStatistics::jump_halos_send_++;
      nid_t const parent_id = ParentIdOfNode(id);

      // buffer size is dependent on jump location
      if (LTI(location) <= LTI(BoundaryLocation::Bottom)) {
        // Plane
        jump_send_counter_plane += UpdateMaterialJumpMpiSend(
            parent_id, requests, id,
            &jump_buffer_plane[jump_send_counter_plane], location, field_type);
      } else if (LTI(location) <= LTI(BoundaryLocation::SouthWest)) {
        // Stick
        jump_send_counter_stick += UpdateMaterialJumpMpiSend(
            parent_id, requests, id,
            &jump_buffer_stick[jump_send_counter_stick], location, field_type);
      } else {
        // Cube
        jump_send_counter_cube += UpdateMaterialJumpMpiSend(
            parent_id, requests, id, &jump_buffer_cube[jump_send_counter_cube],
            location, field_type);
      }
    }
#endif
    }
  }
  if (jump_send_counter_plane > jump_buffer_plane.size()) {
    throw std::out_of_range("jump_buffer_plane");
  }
  if (jump_send_counter_stick > jump_buffer_stick.size()) {
    throw std::out_of_range("jump_buffer_stick");
  }
  if (jump_send_counter_cube > jump_buffer_cube.size()) {
    throw std::out_of_range("jump_buffer_cube");
  }
}

/**
 * @brief Method used to execute the Halo Update for all no-jump boundaries with
 * MPI Communication.
 * @param requests vector of communication request (handle), new handle will be
 * added at the end of the vector.
 * @param boundaries Description of internal boundaries without jump.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void InternalHaloManager::MpiMaterialHaloUpdateNoJump(
    std::vector<MPI_Request> &requests,
    std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const
        &boundaries,
    MaterialFieldType const field_type) {
  for (auto const &boundary : boundaries) {
    nid_t id = std::get<0>(boundary);
    BoundaryLocation location = std::get<1>(boundary);

    switch (std::get<2>(boundary)) {
    case InternalBoundaryType::NoJumpBoundaryMpiSend: {
      CommunicationStatistics::no_jump_halos_send_++;
      UpdateMaterialHaloCellsMpiSend(id, requests, location, field_type);
    } break;
#ifndef PERFORMANCE
    case InternalBoundaryType::NoJumpBoundaryMpiRecv: {
      CommunicationStatistics::no_jump_halos_recv_++;
      UpdateMaterialHaloCellsMpiRecv(id, requests, location, field_type);
    } break;
    default:
      throw std::logic_error(
          "Halo update: BoundaryType does not need MPI, is jump or is unknown");
#else
    default: /* InternalBoundaryType::NoJumpBoundaryMpiRecv */ {
      CommunicationStatistics::no_jump_halos_recv_++;
      UpdateMaterialHaloCellsMpiRecv(id, requests, location, field_type);
    }
#endif
    }
  }
}

/**
 * @brief Method used to execute the Halo Update for all boundaries without MPI
 * Communication.
 * @param boundaries description of internal boundaries, either jump or no jump.
 * @param field_type The decider whether a halo update for conservatives or for
 * prime states is done.
 */
void InternalHaloManager::NoMpiMaterialHaloUpdate(
    std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const
        &boundaries,
    MaterialFieldType const field_type) {
  for (auto const &boundary : boundaries) {
    InternalBoundaryType const type = std::get<2>(boundary);
    nid_t const id = std::get<0>(boundary);
    BoundaryLocation const location = std::get<1>(boundary);

    switch (type) {
    case InternalBoundaryType::JumpBoundaryLocal: {
      UpdateMaterialJumpNoMpi(id, location, field_type);
    } break;
#ifndef PERFORMANCE
    case InternalBoundaryType::NoJumpBoundaryLocal: {
      UpdateMaterialHaloCellsNoMpi(id, location, field_type);
    } break;
    default:
      throw std::logic_error(" Material halo update: unknown boundary type");
#else
    default: /* InternalBoundaryType::NoJumpBoundaryLocal */
      UpdateMaterialHaloCellsNoMpi(id, location, field_type);
#endif
    }
  }
}

/**
 * @brief Executes the Interface tag halo update for all boundaries that do not
 * need communication.
 * @param boundaries Boundaries to be updated.
 * @param buffer_type Level-set buffer type on which the update is done.
 */
void InternalHaloManager::NoMpiInterfaceTagHaloUpdate(
    std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const
        &boundaries,
    InterfaceDescriptionBufferType const buffer_type) {
  for (auto const &boundary : boundaries) {
    nid_t const id = std::get<0>(boundary);
    BoundaryLocation const location = std::get<1>(boundary);
    InternalBoundaryType const type = std::get<2>(boundary);

    if (type != InternalBoundaryType::JumpBoundaryMpiSend) {
      // jump/no_jump is dealt within the boundary condition, but there is no
      // need to project jumps from the parent, therefore it is ignored here.
      UpdateInterfaceTagHaloCellsNoMpi(id, buffer_type, location);
    }
  }
}

/**
 * @brief Executes the interface tag halo update for all boundaries that do need
 * communication.
 * @param boundaries Boundaries to be updated.
 * @param type Level-set buffer type on which the update is done.
 */
void InternalHaloManager::MpiInterfaceTagHaloUpdate(
    std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const
        &boundaries,
    InterfaceDescriptionBufferType const type,
    std::vector<MPI_Request> &requests) {
  for (auto const &boundary : boundaries) {
    BoundaryLocation const location = std::get<1>(boundary);
    nid_t const id = std::get<0>(boundary);

    switch (std::get<2>(boundary)) {
    case InternalBoundaryType::NoJumpBoundaryMpiSend:
      UpdateInterfaceTagHaloCellsMpiSend(id, requests, type, location);
      break;
#ifndef PERFORMANCE
    case InternalBoundaryType::NoJumpBoundaryMpiRecv:
      UpdateInterfaceTagHaloCellsMpiRecv(id, requests, type, location);
      break;
    default:
      throw std::logic_error("BoundaryManager::MpiInterfaceTagHaloUpdate: "
                             "BoundaryType not supported");
#else
    default: /* InternalBoundaryType::NoJumpBoundaryMpiRecv */ {
      UpdateInterfaceTagHaloCellsMpiRecv(id, requests, type, location);
    }
#endif
    }
  }
}

/**
 * @brief Executes the InterfaceHaloUpdate for all boundaries that do not need
 * communication.
 * @param boundaries Boundaries to be updated.
 * @param buffer_type Type of buffer on the levelset block that should be
 * updated.
 */
void InternalHaloManager::NoMpiInterfaceHaloUpdate(
    std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const
        &boundaries,
    InterfaceBlockBufferType const buffer_type) {
  for (auto const &boundary : boundaries) {
    nid_t const id = std::get<0>(boundary);
    BoundaryLocation const location = std::get<1>(boundary);
    InternalBoundaryType const type = std::get<2>(boundary);

    if (type != InternalBoundaryType::JumpBoundaryMpiSend) {
      // jump/no_jump is dealt within the boundary condition, but there is no
      // need to project jumps from the parent, therefore it is ignored here.
      UpdateInterfaceHaloCellsNoMpi(id, buffer_type, location);
    }
  }
}

/**
 * @brief Executes the InterfaceHaloUpdate for all boundaries that do need
 * communication.
 * @param boundaries Boundaries to be updated.
 * @param buffer_type Type of buffer on the levelset block that should be
 * updated.
 */
void InternalHaloManager::MpiInterfaceHaloUpdate(
    std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const
        &boundaries,
    InterfaceBlockBufferType const buffer_type,
    std::vector<MPI_Request> &requests) {
  for (auto const &boundary : boundaries) {
    BoundaryLocation const location = std::get<1>(boundary);
    nid_t const id = std::get<0>(boundary);

    switch (std::get<2>(boundary)) {
    case InternalBoundaryType::NoJumpBoundaryMpiSend: {
      UpdateInterfaceHaloCellsMpiSend(id, requests, buffer_type, location);
    } break;

#ifndef PERFORMANCE
    case InternalBoundaryType::NoJumpBoundaryMpiRecv: {
      UpdateInterfaceHaloCellsMpiRecv(id, requests, buffer_type, location);
    } break;
    default:
      throw std::logic_error("BoundaryManager::MpiInterfaceHaloUpdate: "
                             "BoundaryType not supported");
#else
    default: /* InternalBoundaryType::NoJumpBoundaryLocal */
      UpdateInterfaceHaloCellsMpiRecv(id, requests, buffer_type, location);
#endif
    }
  }
}
