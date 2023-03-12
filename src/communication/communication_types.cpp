//===---------------------- communication_types.cpp -----------------------===//
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
#include "communication_types.h"

#include "block_definitions/field_material_definitions.h"
#include "boundary_condition/boundary_specifications.h"
#include <stdexcept>

/**
 * @brief Constructor creating and allocating the derived MPI_Datatypes.
 */
CommunicationTypes::CommunicationTypes() { CreateDataTypes(); }

/**
 * @brief Destructor releases the allocated Datatypes.
 */
CommunicationTypes::~CommunicationTypes() { FreeTypes(); }

/**
 * @brief Datatypes for 3D Simulations. See Meta function.
 */
void CommunicationTypes::CreateDataTypes() {

  // creates all Datatypes for Halo Updates
  for (unsigned int type = 0; type < 2; type++) {
    // material data types
    for (BoundaryLocation const location : CC::HBS()) {
      // send
      MPI_Type_create_subarray(3, block_size_.data(),
                               halo_size_[LTI(location)].data(),
                               start_indices_halo_send_[LTI(location)].data(),
                               MPI_ORDER_C, type == 0 ? MPI_DOUBLE : MPI_INT8_T,
                               &send_types_[type][LTI(location)]);
      MPI_Type_commit(&send_types_[type][LTI(location)]);

      // recv
      MPI_Type_create_subarray(3, block_size_.data(),
                               halo_size_[LTI(location)].data(),
                               start_indices_halo_recv_[LTI(location)].data(),
                               MPI_ORDER_C, type == 0 ? MPI_DOUBLE : MPI_INT8_T,
                               &recv_types_[type][LTI(location)]);
      MPI_Type_commit(&recv_types_[type][LTI(location)]);
    }

    // ProjectLevel representation of child-memory
    for (unsigned int child = 0; child < CC::NOC(); ++child) {
      MPI_Type_create_subarray(3, block_size_.data(), child_size_.data(),
                               start_index_child_[child].data(), MPI_ORDER_C,
                               type == 0 ? MPI_DOUBLE : MPI_INT8_T,
                               &averaging_send_[type][child]);
      MPI_Type_commit(&averaging_send_[type][child]);
    }
  }

  int const jump_send_start_index[3] = {0, 0, 0};
  MPI_Type_create_subarray(3, halo_size_[LTI(BoundaryLocation::East)].data(),
                           halo_size_[LTI(BoundaryLocation::East)].data(),
                           jump_send_start_index, MPI_ORDER_C, MPI_DOUBLE,
                           &jump_plane_ew_);
  MPI_Type_create_subarray(3, halo_size_[LTI(BoundaryLocation::North)].data(),
                           halo_size_[LTI(BoundaryLocation::North)].data(),
                           jump_send_start_index, MPI_ORDER_C, MPI_DOUBLE,
                           &jump_plane_ns_);
  MPI_Type_create_subarray(3, halo_size_[LTI(BoundaryLocation::Top)].data(),
                           halo_size_[LTI(BoundaryLocation::Top)].data(),
                           jump_send_start_index, MPI_ORDER_C, MPI_DOUBLE,
                           &jump_plane_tb_);

  MPI_Type_commit(&jump_plane_ew_);
  MPI_Type_commit(&jump_plane_ns_);
  MPI_Type_commit(&jump_plane_tb_);

  MPI_Type_create_subarray(
      3, halo_size_[LTI(BoundaryLocation::BottomSouth)].data(),
      halo_size_[LTI(BoundaryLocation::BottomSouth)].data(),
      jump_send_start_index, MPI_ORDER_C, MPI_DOUBLE, &jump_stick_x_);
  MPI_Type_create_subarray(
      3, halo_size_[LTI(BoundaryLocation::BottomWest)].data(),
      halo_size_[LTI(BoundaryLocation::BottomWest)].data(),
      jump_send_start_index, MPI_ORDER_C, MPI_DOUBLE, &jump_stick_y_);
  MPI_Type_create_subarray(
      3, halo_size_[LTI(BoundaryLocation::SouthWest)].data(),
      halo_size_[LTI(BoundaryLocation::SouthWest)].data(),
      jump_send_start_index, MPI_ORDER_C, MPI_DOUBLE, &jump_stick_z_);

  MPI_Type_commit(&jump_stick_x_);
  MPI_Type_commit(&jump_stick_y_);
  MPI_Type_commit(&jump_stick_z_);

  MPI_Type_create_subarray(
      3, halo_size_[LTI(BoundaryLocation::EastNorthTop)].data(),
      halo_size_[LTI(BoundaryLocation::EastNorthTop)].data(),
      jump_send_start_index, MPI_ORDER_C, MPI_DOUBLE, &jump_cube_);
  MPI_Type_commit(&jump_cube_);

  int const tc_per_conservative = CC::TCX() * CC::TCY() * CC::TCZ();
  MPI_Type_contiguous(tc_per_conservative, MPI_DOUBLE, &single_conservatives_);
  int const tc_per_jump = MF::ANOE() * CC::ICY() * CC::ICZ();
  MPI_Type_contiguous(tc_per_jump, MPI_DOUBLE, &single_boundary_jump_);

  MPI_Type_commit(&single_conservatives_);
  MPI_Type_commit(&single_boundary_jump_);
}

/**
 * @brief Frees MPI datatypes used in three dimensional simulations. See Meta
 * function.
 */
void CommunicationTypes::FreeTypes() {

  for (unsigned int type = 0; type < 2; type++) {
    // material data types
    for (BoundaryLocation location : CC::HBS()) {
      MPI_Type_free(&send_types_[type][LTI(location)]);
      MPI_Type_free(&recv_types_[type][LTI(location)]);
    }

    // ProjectLevel representation of child-memory
    for (unsigned int i = 0; i < CC::NOC(); ++i) {
      MPI_Type_free(&averaging_send_[type][i]);
    }
  }

  // Jump Halo data types
  MPI_Type_free(&jump_plane_ew_);
  MPI_Type_free(&jump_plane_ns_);
  MPI_Type_free(&jump_plane_tb_);
  MPI_Type_free(&jump_stick_x_);
  MPI_Type_free(&jump_stick_y_);
  MPI_Type_free(&jump_stick_z_);
  MPI_Type_free(&jump_cube_);

  // Whole block-struct
  MPI_Type_free(&single_conservatives_);
  MPI_Type_free(&single_boundary_jump_);
}

/**
 * @brief Gives the MPI Datatype to send data into in an internal no-jump
 * exchange.
 * @param location The direction of the halo to be sent into.
 * @param datatype Mpi Datatype that should be received.
 * @return The Datatype needed to send into the right sub_array.
 */
MPI_Datatype
CommunicationTypes::RecvDatatype(BoundaryLocation const location,
                                 DatatypeForMpi const datatype) const {
  return recv_types_[DTI(datatype)][LTI(location)];
}

/**
 * @brief Gives the MPI Datatype to send data from in an internal cell.
 * @param location The direction of the halo to be sent from.
 * @param datatype Mpi Datatype that should be sent.
 * @return The Datatype needed to send from the right sub_array.
 */
MPI_Datatype
CommunicationTypes::SendDatatype(BoundaryLocation const location,
                                 DatatypeForMpi const datatype) const {
  return send_types_[DTI(datatype)][LTI(location)];
}

/**
 * @brief Gives a MPI Datatype to send the full set of conservatives.
 * @return The datatype for conservative buffer communication.
 */
MPI_Datatype CommunicationTypes::ConservativesDatatype() const {
  return single_conservatives_;
}

/**
 * @brief Gives a MPI Datatype to send the full set of jump buffers.
 * @return The datatype for jump surface buffer communication.
 */
MPI_Datatype CommunicationTypes::JumpSurfaceDatatype() const {
  return single_boundary_jump_;
}

/**
 * @brief Gives the MPI Datatype to send the tags from a child into a parent.
 * @param child_position The position of the child among the siblings. See Node
 * class for more details.
 * @param datatype Mpi Datatype that should be sent.
 * @return The subarray into which the data is to be averaged.
 */
MPI_Datatype
CommunicationTypes::AveragingSendDatatype(unsigned int const child_position,
                                          DatatypeForMpi const datatype) const {
  if (child_position < CC::NOC()) {
    return averaging_send_[DTI(datatype)][child_position];
  }
  throw std::logic_error("Child position not possible: " +
                         std::to_string(child_position));
}

/**
 * @brief Gives the correct datatype for a jump boundary at the given location.
 * @param location The location identifier.
 * @return The correct MPI_Datatype.
 */
MPI_Datatype CommunicationTypes::JumpPlaneSendDatatype(
    BoundaryLocation const location) const {

  switch (location) {
  case BoundaryLocation::East:
  case BoundaryLocation::West:
    return jump_plane_ew_;
  case BoundaryLocation::North:
  case BoundaryLocation::South:
    return jump_plane_ns_;
  case BoundaryLocation::Top:
  case BoundaryLocation::Bottom:
    return jump_plane_tb_;
  case BoundaryLocation::BottomNorth:
  case BoundaryLocation::BottomSouth:
  case BoundaryLocation::TopNorth:
  case BoundaryLocation::TopSouth:
    return jump_stick_x_;
  case BoundaryLocation::BottomEast:
  case BoundaryLocation::BottomWest:
  case BoundaryLocation::TopEast:
  case BoundaryLocation::TopWest:
    return jump_stick_y_;
  case BoundaryLocation::NorthEast:
  case BoundaryLocation::NorthWest:
  case BoundaryLocation::SouthEast:
  case BoundaryLocation::SouthWest:
    return jump_stick_z_;
  default:
    return jump_cube_;
  }
}
