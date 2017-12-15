/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*  \\\\                                                                                  *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA                                                                                 *
* Copyright (c) 2017 Nikolaus A. Adams and contributors (see AUTHORS list)               *
* All rights reserved.                                                                   *
*                                                                                        *
* Chair of Aerodynamics and Fluid Mechanics                                              *
* Technical University of Munich                                                         *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
* This project has received funding from the European Reseach Council (ERC)              *
* under the European Union's Horizon 2020 research and innovation programme              *
* (grant agreement No 667483).                                                           *
*                                                                                        *
* ERC Advanced Grant No 667483, Prof. Dr. Nikolaus A. Adams:                             *
* "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"      *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Redistribution and use in source and binary forms, with or without                     *
* modification, are permitted provided that the following conditions are met:            *
*                                                                                        *
* 1. Redistributions of source code must retain the above copyright notice,              *
*    this list of conditions and the following disclaimer.                               *
*                                                                                        *
* 2. Redistributions in binary form must reproduce the above copyright notice            *
*    this list of conditions and the following disclaimer in the documentation           *
*    and/or other materials provided with the distribution.                              *
*                                                                                        *
* 3. Neither the name of the copyright holder nor the names of its                       *
*    contributors may be used to endorse or promote products derived from this           *
*    software without specific prior written permission.                                 *
*                                                                                        *
* 4. Any redistribution of substantial fractions of the code as a                        *
*    different project should preserve the word ALPACA in the name                       *
*    of the code                                                                         *
*                                                                                        *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"            *
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE              *
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE            *
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE              *
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                    *
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF                   *
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS               *
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN                *
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)                *
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE            *
* POSSIBILITY OF SUCH DAMAGE.                                                            *
*                                                                                        *
* Please note, several third-party tools are used within the ALPACA code under           *
* their own license agreement.                                                           *
*                                                                                        *
* 1. xdmf_writer        : Licensed by Technische Universitaet Muenchen                   *
*                         See 'COPYING_XDMF_WRITER' for more information.                *
*                                                                                        *
* 2. tiny_xml           : This software is provided 'as-is', without any express or      *
*                         implied warranty. In no event will the authors be held         *
*                         liable for any damages arising from the use of this software.  *
*                         See COPYING_TINY_XMLfor more information.                      *
*                                                                                        *
* 3. expression_toolkit : Free use of The C++ Mathematical Expression Toolkit Library is *
*                         permitted under the guidelines and in accordance with the most *
*                         current version of the Common Public License.                  *
*                         http://www.opensource.org/licenses/cpl1.0.php                  *
*                         See COPYING_EXPRESSION_TOOLKITfor more information.            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* AUTHORS                                                                                *
*                                                                                        *
*   Prof. Dr. Nikolaus A. Adams                                                          *
*                                                                                        *
*   Dr. Stefan Adami                                                                     *
*   Vladimir Bogdanov                                                                    *
*   Nico Fleischmann                                                                     *
*   Nils Hoppe                                                                           *
*   Naeimeh Hosseini                                                                     *
*   Jakob Kaiser                                                                         *
*   Aleksandr Lunkov                                                                     *
*   Thomas Paula                                                                         *
*   Josef Winter                                                                         *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
*   nanoshock@aer.mw.tum.de                                                              *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, December 15th 2017                                                             *
*                                                                                        *
*****************************************************************************************/

#include "internal_boundary_condition.h"

#include <memory>
#include <stdexcept>
#include "node.h"
#include "multi_resolution/multi_resolution.h"

/**
 * @brief See base class constructor.
 */
InternalBoundaryCondition::InternalBoundaryCondition(const BoundaryLocation location, const Tree& tree, const TopologyManager& topology, const std::uint64_t host_id) :
    BoundaryCondition(location, tree, topology, host_id),
    neigbor_id_(Node::GetNeighborId(host_id_,location_))
{
}

/**
 * @brief Updates halos without neighbor on the same level. In order to fill the halo cells a prediction of the parent values is
 *        performed. All communication is blocking.
 */
void InternalBoundaryCondition::UpdateJump() {

  Block& host_block = tree_.GetNodeWithId(host_id_)->GetBlock();

  unsigned int x_start_index = 0;
  unsigned int y_start_index = 0;
  unsigned int z_start_index = 0;
  unsigned int x_end_index = CC::TCX();
  unsigned int y_end_index = CC::TCY();
  unsigned int z_end_index = CC::TCZ();

  switch (location_) {
    case BoundaryLocation::eEast : {
       x_start_index = CC::FHH();
    }
    break;
    case BoundaryLocation::eWest : {
        x_end_index = CC::HS();
    }
    break;
    case BoundaryLocation::eNorth : {
        y_start_index = CC::FHH();
    }
    break;
    case BoundaryLocation::eSouth : {
        y_end_index = CC::HS();
    }
    break;
    case BoundaryLocation::eTop : {
        z_start_index = CC::FHH();
    }
    break;
    case BoundaryLocation::eBottom : {
        z_end_index = CC::HS();
    }
    break;
    default:
        throw std::invalid_argument("Error in Internal BC - Location not found");
    break;
  }

  if(topology_.NodeIsOnMyRank(Node::ParentIdOfNode(host_id_))) {
    std::shared_ptr<Node> parent = tree_.GetNodeWithId(Node::ParentIdOfNode(host_id_));
    Block& parent_block = parent->GetBlock();

    MultiResolution::Prediction(parent_block.GetRightHandSideBuffer(0), host_block.GetRightHandSideBuffer(0), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);
    MultiResolution::Prediction(parent_block.GetRightHandSideBuffer(1), host_block.GetRightHandSideBuffer(1), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);
    MultiResolution::Prediction(parent_block.GetRightHandSideBuffer(2), host_block.GetRightHandSideBuffer(2), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);
    MultiResolution::Prediction(parent_block.GetRightHandSideBuffer(3), host_block.GetRightHandSideBuffer(3), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);
    MultiResolution::Prediction(parent_block.GetRightHandSideBuffer(4), host_block.GetRightHandSideBuffer(4), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);

  } else {
    // MPI
      MPI_Datatype recv_type = topology_.GetDatatypeForJumpBoundary(std::make_pair(host_id_,location_));

      int sender_rank = topology_.GetRankOfNode(Node::ParentIdOfNode(host_id_));

             double parent_rho[CC::TCX()][CC::TCY()][CC::TCZ()];
          double parent_energy[CC::TCX()][CC::TCY()][CC::TCZ()];
      double parent_x_momentum[CC::TCX()][CC::TCY()][CC::TCZ()];
      double parent_y_momentum[CC::TCX()][CC::TCY()][CC::TCZ()];
      double parent_z_momentum[CC::TCX()][CC::TCY()][CC::TCZ()];

      MPI_Recv(&parent_rho,1,recv_type,sender_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&parent_energy,1,recv_type,sender_rank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&parent_x_momentum,1,recv_type,sender_rank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&parent_y_momentum,1,recv_type,sender_rank,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&parent_z_momentum,1,recv_type,sender_rank,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

      MultiResolution::Prediction(parent_rho,        host_block.GetRightHandSideBuffer(0), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);
      MultiResolution::Prediction(parent_energy,     host_block.GetRightHandSideBuffer(1), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);
      MultiResolution::Prediction(parent_x_momentum, host_block.GetRightHandSideBuffer(2), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);
      MultiResolution::Prediction(parent_y_momentum, host_block.GetRightHandSideBuffer(3), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);
      MultiResolution::Prediction(parent_z_momentum, host_block.GetRightHandSideBuffer(4), host_id_, x_start_index, x_end_index, y_start_index, y_end_index, z_start_index, z_end_index);
  }
}

/**
 * @brief Fill the halo cells with the facing domain values of the neighbor block.
 * @param requests All communication is non blocking, therefore requests need to be propagted to this outside container parameter
 */
void InternalBoundaryCondition::UpdateNoJump(std::vector<MPI_Request>& requests) {

 Block& host_block = tree_.GetNodeWithId(host_id_)->GetBlock(); // <- Dirty Fix, should be handled when proper Boundaries are inserted

 int rank_of_neighbor = topology_.GetRankOfNode(neigbor_id_);

  if(topology_.GetRankOfNode(host_id_) != rank_of_neighbor) {
    //Do MPI
    int recv_base_tag = topology_.GetTagForNode(host_id_) << 4;
    int send_base_tag = topology_.GetTagForNode(neigbor_id_) << 4;

    MPI_Datatype send_type;
    MPI_Datatype recv_type;

    send_type = topology_.NoJumpSliceDataType(location_);
    recv_type = topology_.NoJumpBoundaryDataType(location_);

    // Avoids sending from wrong swide if block is sourrounded by same rank, i.e. #1 - #2 - #1
    switch (location_) {
      case BoundaryLocation::eEast : {
        send_base_tag += 8;
      }
      break;
      case BoundaryLocation::eWest : {
        recv_base_tag += 8;
      }
      break;
      case BoundaryLocation::eNorth : {
        send_base_tag += 8;
      }
      break;
      case BoundaryLocation::eSouth : {
        recv_base_tag += 8;
      }
      break;
      case BoundaryLocation::eTop : {
        send_base_tag += 8;
      }
      break;
      case BoundaryLocation::eBottom : {
        recv_base_tag += 8;
      }
      break;
      default:
        throw std::invalid_argument("Error in Internal BC - Location not found");
      break;
      }

      for(unsigned int e=0;e<CC::NoEq();e++){
          //Due to MPI's C Bindings we need to to this a little cumbersome
          requests.push_back(MPI_Request());
          MPI_Isend(host_block.GetRightHandSideBuffer(e),1,send_type, rank_of_neighbor, send_base_tag+e,MPI_COMM_WORLD,&requests.back());
          requests.push_back(MPI_Request());
          MPI_Irecv(host_block.GetRightHandSideBuffer(e),1,recv_type, rank_of_neighbor, recv_base_tag+e,MPI_COMM_WORLD,&requests.back()); // Caution uses "back()" again
      }
  } else {
    //Do non MPI
    Block& partner_block = tree_.GetNodeWithId(neigbor_id_)->GetBlock();
    unsigned int offset = CC::FHH() - CC::HS();
    switch (location_) {
      case BoundaryLocation::eEast : {
        for(unsigned int e=0; e < CC::NoEq(); ++e) {
            double (&host_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            double (&partner_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = partner_block.GetRightHandSideBuffer(e);
            for(unsigned int i = CC::FHH(); i < CC::TCX(); ++i ) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        host_cells[i][j][k] = partner_cells[i-offset][j][k];
                    }
                }
            }
        }
      }
      break;
      case BoundaryLocation::eWest : {
        for(unsigned int e=0; e < CC::NoEq(); ++e) {
            double (&host_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            double (&partner_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = partner_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::HS(); ++i ) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        host_cells[i][j][k] = partner_cells[i+offset][j][k];
                    }
                }
            }
        }
      }
      break;
      case BoundaryLocation::eNorth : {
        for(unsigned int e=0; e < CC::NoEq(); ++e) {
            double (&host_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            double (&partner_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = partner_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::TCX(); ++i ) {
                for(unsigned int j = CC::FHH(); j < CC::TCY(); ++j ) {
                    for(unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        host_cells[i][j][k] = partner_cells[i][j-offset][k];
                    }
                }
            }
        }
      }
      break;
      case BoundaryLocation::eSouth : {
        for(unsigned int e=0; e < CC::NoEq(); ++e) {
            double (&host_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            double (&partner_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = partner_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::TCX(); ++i ) {
                for(unsigned int j = 0; j < CC::HS(); ++j ) {
                    for(unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        host_cells[i][j][k] = partner_cells[i][j+offset][k];
                    }
                }
            }
        }
      }
      break;
      case BoundaryLocation::eTop : {
        for(unsigned int e=0; e < CC::NoEq(); ++e) {
            double (&host_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            double (&partner_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = partner_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::TCX(); ++i ) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int k = CC::FHH(); k < CC::TCZ(); ++k ) {
                        host_cells[i][j][k] = partner_cells[i][j][k-offset];
                    }
                }
            }
        }
      }
      break;
      case BoundaryLocation::eBottom : {
        for(unsigned int e=0; e < CC::NoEq(); ++e) {
            double (&host_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            double (&partner_cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = partner_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::TCX(); ++i ) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int k = 0; k < CC::HS(); ++k ) {
                        host_cells[i][j][k] = partner_cells[i][j][k+offset];
                    }
                }
            }
        }
      }
      break;
      default:
        throw std::invalid_argument("Error in Internal BC - Location not found");
      break;
    }
  }

}
/**
 * @brief See base class.
 * @param requests Outside Container for for non-blocking MPI communication requests. Only used in no-jump case, however, the state is not known from outside
 *        therefore the parameter is needed.
 */
void InternalBoundaryCondition::UpdateHaloCells(std::vector<MPI_Request>& requests) {
  if(topology_.NodeExists(neigbor_id_)) {
    UpdateNoJump(requests);
  } else {
    UpdateJump();
  }
}

/**
 * @brief See base class
 * @return eInternal, always.
 */
BoundaryType InternalBoundaryCondition::GetType() const {
    return BoundaryType::eInternal;
}
