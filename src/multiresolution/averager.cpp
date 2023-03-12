//===-------------------------- averager.cpp ------------------------------===//
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
#include "multiresolution/averager.h"
#include "communication/communication_statistics.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "multiresolution/multiresolution.h"
#include "topology/id_information.h"
#include "user_specifications/debug_and_profile_setup.h"

namespace {
/**
 * @brief Filters zero entries from the input.
 * @param descending_values Must be unique and ordered descending.
 * @return copy of the input vector with the zero entries cut.
 */
std::vector<unsigned int> DescendingVectorWithoutZero(
    std::vector<unsigned int> const &descending_values) {
  return std::vector<unsigned int>(std::cbegin(descending_values),
                                   std::find(std::cbegin(descending_values),
                                             std::cend(descending_values), 0));
}
} // namespace

/**
 * @brief Default constructor.
 * @param topology The underlying topology.
 * @param communicator The communcation manger used for sending and receiving
 * remote data.
 * @param tree The tree to be worked on.
 */
Averager::Averager(TopologyManager const &topology,
                   CommunicationManager &communicator, Tree &tree)
    : topology_(topology), communicator_(communicator), tree_(tree) {
  // Empty besides initializer list.
}

/**
 * @brief Fill parents of nodes on the given level via conservative average
 * operations.
 * @param child_levels_descending The levels of the children holding the data to
 * be averaged.
 */
void Averager::AverageMaterial(
    std::vector<unsigned int> const &child_levels_descending) const {

  for (unsigned int const child_level :
       DescendingVectorWithoutZero(child_levels_descending)) {

    std::vector<std::tuple<nid_t, int, int>>
        child_parent_rank_relations; // id, rank-child, rank-parent
    std::vector<nid_t> no_mpi_list;
    unsigned int send_counter = 0;

    // sort by mpi necessary or not
    for (nid_t const child_id : topology_.IdsOnLevel(child_level)) {
      nid_t const parent_id = ParentIdOfNode(child_id);
      int const rank_of_child = topology_.GetRankOfNode(child_id);
      int const rank_of_parent = topology_.GetRankOfNode(parent_id);

      if (rank_of_child == communicator_.MyRankId() &&
          rank_of_parent == communicator_.MyRankId()) {
        no_mpi_list.push_back(child_id);
      } else if (rank_of_child == communicator_.MyRankId() ||
                 rank_of_parent == communicator_.MyRankId()) {
        child_parent_rank_relations.push_back(
            std::make_tuple(child_id, rank_of_child, rank_of_parent));
        if (rank_of_child == communicator_.MyRankId()) {
          send_counter += topology_.GetMaterialsOfNode(child_id)
                              .size(); // needed for buffer size
        }
      }
    }

    std::vector<Conservatives> send_buffer_parent(send_counter);
    send_counter = 0; // is reused
    std::vector<MPI_Request> requests;

    for (auto const &[child_id, rank_of_child, rank_of_parent] :
         child_parent_rank_relations) {
      if (rank_of_child == communicator_.MyRankId() &&
          rank_of_parent != communicator_.MyRankId()) {
        Node const &child = tree_.GetNodeWithId(child_id);
        unsigned int const pos = PositionOfNodeAmongSiblings(child_id);
        for (auto const material : topology_.GetMaterialsOfNode(child_id)) {
          Multiresolution::Average(
              child.GetPhaseByMaterial(material).GetRightHandSideBuffer(),
              send_buffer_parent.at(send_counter), child_id);
          communicator_.Send(
              &send_buffer_parent.at(send_counter), MF::ANOE(),
              communicator_.AveragingSendDatatype(pos, DatatypeForMpi::Double),
              rank_of_parent, requests);
          send_counter++;
          if constexpr (DP::Profile()) {
            CommunicationStatistics::average_level_send_++;
          }
        }
      } else if (rank_of_child != communicator_.MyRankId() &&
                 rank_of_parent == communicator_.MyRankId()) {
        nid_t const parent_id = ParentIdOfNode(child_id);
        Node &parent = tree_.GetNodeWithId(parent_id);
        unsigned int const pos = PositionOfNodeAmongSiblings(child_id);
        for (auto const material : topology_.GetMaterialsOfNode(child_id)) {
          communicator_.Recv(
              &parent.GetPhaseByMaterial(material).GetRightHandSideBuffer(),
              MF::ANOE(),
              communicator_.AveragingSendDatatype(pos, DatatypeForMpi::Double),
              rank_of_child, requests);
          if constexpr (DP::Profile()) {
            CommunicationStatistics::average_level_recv_++;
          }
        }
      }
    }

    for (nid_t const child_id : no_mpi_list) {
      nid_t const parent_id = ParentIdOfNode(child_id);
      Node &parent = tree_.GetNodeWithId(parent_id);
      Node const &child = tree_.GetNodeWithId(child_id);
      for (auto const material : topology_.GetMaterialsOfNode(child_id)) {
        Multiresolution::Average(
            child.GetPhaseByMaterial(material).GetRightHandSideBuffer(),
            parent.GetPhaseByMaterial(material).GetRightHandSideBuffer(),
            child_id);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
  } // child_level
}

/**
 * @brief Projects the interface tag down to the parent levels.
 * @param levels_with_updated_parents_descending The child levels whose parents
 * were updated in descending order.
 */
void Averager::AverageInterfaceTags(
    std::vector<unsigned int> const &levels_with_updated_parents_descending)
    const {

  for (unsigned int const child_level :
       DescendingVectorWithoutZero(levels_with_updated_parents_descending)) {
    std::vector<std::tuple<nid_t, int, int>>
        child_parent_rank_relations; // child-id, rank-child, rank-parent
    std::vector<nid_t> no_mpi_same_rank;
    std::vector<nid_t> no_mpi_uniform_child;
    unsigned int send_counter = 0;

    // sort into list
    for (nid_t const child_id : topology_.IdsOnLevel(child_level)) {
      nid_t const parent_id = ParentIdOfNode(child_id);

      // if the parent is NOT multi, the child is neither and we do not have to
      // do anything with the tags
      if (topology_.IsNodeMultiPhase(parent_id)) {
        int const rank_of_child = topology_.GetRankOfNode(child_id);
        int const rank_of_parent = topology_.GetRankOfNode(parent_id);
        if (rank_of_child == communicator_.MyRankId() &&
            rank_of_parent == communicator_.MyRankId()) {
          // Non MPI Averaging
          no_mpi_same_rank.emplace_back(child_id);

        } else if (rank_of_child == communicator_.MyRankId() &&
                   rank_of_parent != communicator_.MyRankId()) {
          // if the child is not multi, the parent figures its uniform tag out
          // by its own
          if (topology_.IsNodeMultiPhase(child_id)) {
            child_parent_rank_relations.push_back(
                std::make_tuple(child_id, rank_of_child, rank_of_parent));
            send_counter++;
          }
        } else if (rank_of_child != communicator_.MyRankId() &&
                   rank_of_parent == communicator_.MyRankId()) {
          if (topology_.IsNodeMultiPhase(child_id)) {
            child_parent_rank_relations.push_back(
                std::make_tuple(child_id, rank_of_child, rank_of_parent));
          } else {
            // parent is updated without communication
            no_mpi_uniform_child.emplace_back(child_id);
          }
        }
      } // if NodeIsMulti
    }   // children on child_level

    std::vector<InterfaceTagBundle> send_buffer_parent(
        send_counter); // buffer necessary for asynchronous send of averaged
                       // values
    std::vector<MPI_Request> requests;
    send_counter = 0;

    // MPI Communications
    for (auto const &[child_id, rank_of_child, rank_of_parent] :
         child_parent_rank_relations) {
      if (rank_of_child == communicator_.MyRankId() &&
          rank_of_parent != communicator_.MyRankId()) {
        Node const &child = tree_.GetNodeWithId(child_id);
        Multiresolution::PropagateCutCellTagsFromChildIntoParent(
            child.GetInterfaceTags<
                InterfaceDescriptionBufferType::Reinitialized>(),
            send_buffer_parent.at(send_counter).interface_tags_, child_id);
        int const pos = PositionOfNodeAmongSiblings(child_id);
        communicator_.Send(
            &send_buffer_parent.at(send_counter), 1,
            communicator_.AveragingSendDatatype(pos, DatatypeForMpi::Byte),
            rank_of_parent, requests);
        send_counter++;
      } else if (rank_of_child != communicator_.MyRankId() &&
                 rank_of_parent == communicator_.MyRankId()) {
        // Recv
        nid_t const parent_id = ParentIdOfNode(child_id);
        Node &parent = tree_.GetNodeWithId(parent_id);
        int const pos = PositionOfNodeAmongSiblings(child_id);
        communicator_.Recv(
            &parent.GetInterfaceTags<
                InterfaceDescriptionBufferType::Reinitialized>(),
            1, communicator_.AveragingSendDatatype(pos, DatatypeForMpi::Byte),
            rank_of_child, requests);
      }
    }

    // figure out the uniform tag of the child without MPI
    for (auto child_id : no_mpi_uniform_child) {
      nid_t const parent_id = ParentIdOfNode(child_id);
      Node &parent = tree_.GetNodeWithId(parent_id);
      std::int8_t const uniform_tag =
          MaterialSignCapsule::SignOfMaterial(
              topology_.GetMaterialsOfNode(child_id).back()) *
          ITTI(IT::BulkPhase);
      Multiresolution::PropagateUniformTagsFromChildIntoParent(
          uniform_tag,
          parent.GetInterfaceTags<
              InterfaceDescriptionBufferType::Reinitialized>(),
          child_id);
    }

    // Non MPI Averaging
    for (auto const child_id : no_mpi_same_rank) {
      nid_t const parent_id = ParentIdOfNode(child_id);
      Node &parent = tree_.GetNodeWithId(parent_id);
      if (topology_.IsNodeMultiPhase(child_id)) {
        Node const &child = tree_.GetNodeWithId(child_id);
        Multiresolution::PropagateCutCellTagsFromChildIntoParent(
            child.GetInterfaceTags<
                InterfaceDescriptionBufferType::Reinitialized>(),
            parent.GetInterfaceTags<
                InterfaceDescriptionBufferType::Reinitialized>(),
            child_id);
      } else {
        std::int8_t const uniform_tag =
            MaterialSignCapsule::SignOfMaterial(
                topology_.GetMaterialsOfNode(child_id).back()) *
            ITTI(IT::BulkPhase);
        Multiresolution::PropagateUniformTagsFromChildIntoParent(
            uniform_tag,
            parent.GetInterfaceTags<
                InterfaceDescriptionBufferType::Reinitialized>(),
            child_id);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
  }
}

/**
 * @brief Fill parents of nodes on the given level via projection operations for
 * all parameters.
 * @param child_levels_descending The levels of the children holding the data to
 * be projected.
 */
void Averager::AverageParameters(
    std::vector<unsigned int> const child_levels_descending) const {
  for (unsigned int child_level : child_levels_descending) {
    // if level 0 should be reached, break the loop as there is no parent level
    // available
    if (child_level == 0)
      break;

    std::vector<MPI_Request> requests;
    std::vector<std::tuple<nid_t, int, int>>
        child_parent_rank_relations; // id, rank-child, rank-parent
    unsigned int sendcounter = 0;
    std::vector<nid_t> no_mpi_list;

    // sort by mpi necessary or not
    for (nid_t const child_id : topology_.IdsOnLevel(child_level)) {
      nid_t const parent_id = ParentIdOfNode(child_id);
      int const rank_of_child = topology_.GetRankOfNode(child_id);
      int const rank_of_parent = topology_.GetRankOfNode(parent_id);

      if (rank_of_child == communicator_.MyRankId() &&
          rank_of_parent == communicator_.MyRankId()) {
        // Non MPI Projection

        no_mpi_list.push_back(child_id);
      } else if (rank_of_child == communicator_.MyRankId() ||
                 rank_of_parent == communicator_.MyRankId()) {
        // MPI Projection

        child_parent_rank_relations.push_back(
            std::make_tuple(child_id, rank_of_child, rank_of_parent));
        if (rank_of_child == communicator_.MyRankId()) {
          sendcounter += topology_.GetMaterialsOfNode(child_id)
                             .size(); // needed for buffer size
        }
      }
    } // nodes on level

    std::vector<Parameters> send_buffer_parent(sendcounter);
    sendcounter = 0;

    for (auto const &rank_relation : child_parent_rank_relations) {
      nid_t const child_id = std::get<0>(rank_relation);
      int const rank_of_child = std::get<1>(rank_relation);
      int const rank_of_parent = std::get<2>(rank_relation);

      if (rank_of_child == communicator_.MyRankId() &&
          rank_of_parent != communicator_.MyRankId()) {
        // MPI_Send
        Node const &child = tree_.GetNodeWithId(child_id);
        unsigned int const pos = PositionOfNodeAmongSiblings(child_id);

        for (auto const material : topology_.GetMaterialsOfNode(child_id)) {
          Multiresolution::Average(
              child.GetPhaseByMaterial(material).GetParameterBuffer(),
              send_buffer_parent.at(sendcounter), child_id);
          communicator_.Send(
              &send_buffer_parent.at(sendcounter), MF::ANOPA(),
              communicator_.AveragingSendDatatype(pos, DatatypeForMpi::Double),
              rank_of_parent, requests);
          sendcounter++;

          if constexpr (DP::Profile()) {
            CommunicationStatistics::average_level_send_++;
          }
        }
      } else if (rank_of_child != communicator_.MyRankId() &&
                 rank_of_parent == communicator_.MyRankId()) {
        // MPI_Recv
        nid_t const parent_id = ParentIdOfNode(child_id);
        Node &parent = tree_.GetNodeWithId(parent_id);
        unsigned int const pos = PositionOfNodeAmongSiblings(child_id);

        for (auto const material : topology_.GetMaterialsOfNode(child_id)) {
          communicator_.Recv(
              &parent.GetPhaseByMaterial(material).GetParameterBuffer(),
              MF::ANOPA(),
              communicator_.AveragingSendDatatype(pos, DatatypeForMpi::Double),
              rank_of_child, requests);
          if constexpr (DP::Profile()) {
            CommunicationStatistics::average_level_recv_++;
          }
        }
      }
    }

    for (nid_t const child_id : no_mpi_list) {
      nid_t const parent_id = ParentIdOfNode(child_id);
      Node &parent = tree_.GetNodeWithId(parent_id);
      Node const &child = tree_.GetNodeWithId(child_id);

      for (auto const material : topology_.GetMaterialsOfNode(child_id)) {
        Multiresolution::Average(
            child.GetPhaseByMaterial(material).GetParameterBuffer(),
            parent.GetPhaseByMaterial(material).GetParameterBuffer(), child_id);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
  } // child_level
}
