//===--------------------- communication_manager.h ------------------------===//
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
#ifndef COMMUNICATION_MANAGER_H
#define COMMUNICATION_MANAGER_H

#include "communication/communication_types.h"
#include "communication/exchange_types.h"
#include "internal_boundary_types.h"
#include "topology/topology_manager.h"
#include "topology/tree.h"
#include <numeric>
#include <vector>

/**
 * @brief The CommunicationManager class provides the functionality for
 * communicating data between nodes and ranks. Furthermore, it holds the
 * neighbor relations between nodes and external relations for external node
 * boundaries.
 */
class CommunicationManager : public CommunicationTypes {
  TopologyManager &topology_;
  unsigned int const maximum_level_;
  int const my_rank_id_;
  int const mpi_tag_ub_;
  std::vector<unsigned int> partner_tag_map_;

  /**
   * Schlachtplan:
   * Make Topology -> unordered_map
   * When remeshing -> go through topology and give two vectors of From-To node
   * pairs and the ranks of each: The first list (send_list) has the local node
   * as from The second list (receive_list) has the local node as to Sort lists
   * first by sender than by receiver (or vice versa) Go through lists and start
   * counting tags (as currently done) for each rank combination. Save the
   * resulting To-From-Tag tuple. Look-up when needed.
   */

  // Cache for Halo Update Pattern
  std::vector<
      std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>>>
      internal_boundaries_;
  std::vector<
      std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>>>
      internal_boundaries_mpi_;
  std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>>
      internal_multi_boundaries_;
  std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>>
      internal_multi_boundaries_mpi_;
  std::vector<
      std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>>>
      internal_boundaries_jump_;
  std::vector<
      std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>>>
      internal_boundaries_jump_mpi_;
  std::vector<std::vector<std::tuple<nid_t, BoundaryLocation>>>
      external_boundaries_;
  std::vector<std::tuple<nid_t, BoundaryLocation>> external_multi_boundaries_;

  // Vector holding flags dor each level that the lists have been created
  // successfully
  std::vector<bool> boundaries_valid_;

  // Three values for three dimensions, even if only one dimension is simulated
  std::vector<std::array<unsigned int, 3>>
      jump_send_count_; // [0]: Plane, [1]: Stick, [2]: cube

  // Function that gives all neighbor-location and external-location relations
  // for a given global node
  void NeighborsOfNode(
      nid_t const global_id,
      std::vector<std::tuple<nid_t, BoundaryLocation>>
          &nodes_internal_boundaries,
      std::vector<std::tuple<nid_t, BoundaryLocation>> &external_boundaries);

public:
  CommunicationManager() = delete;
  explicit CommunicationManager(TopologyManager &topology,
                                unsigned int const maximum_level);
  ~CommunicationManager() = default;
  CommunicationManager(CommunicationManager const &) = delete;
  CommunicationManager &operator=(CommunicationManager const &) = delete;
  CommunicationManager(CommunicationManager &&) = delete;
  CommunicationManager &operator=(CommunicationManager &&) = delete;

  // Function to fill the lists holding the relation to neighbor nodes and
  // external boundaries
  void GenerateNeighborRelationForHaloUpdate(unsigned int const level);

  // return functions for the relation lists
  std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const &
  InternalBoundariesJumpMpi(unsigned level) const;
  std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const &
  InternalBoundariesJump(unsigned level) const;
  std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const &
  InternalBoundariesMpi(unsigned level) const;
  std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const &
  InternalBoundaries(unsigned level) const;
  std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const &
  InternalMultiBoundariesMpi() const;
  std::vector<std::tuple<nid_t, BoundaryLocation, InternalBoundaryType>> const &
  InternalMultiBoundaries() const;
  std::vector<std::tuple<nid_t, BoundaryLocation>> const &
  ExternalBoundaries(unsigned level) const;
  std::vector<std::tuple<nid_t, BoundaryLocation>> const &
  ExternalMultiBoundaries() const;

  // Functions to get the status of the list creations and to empty the flags to
  // regenerate the lists
  bool AreBoundariesValid(unsigned level) const;
  void InvalidateCache();

  // Returns the counter for jump boundaries for the different exchange types
  unsigned int JumpSendCount(unsigned int const level, ExchangeType const type);

  // Send and receive function to buffer data between nodes and ranks
  int Send(void const *buffer, int const count, MPI_Datatype const datatype,
           int const destination_rank, std::vector<MPI_Request> &requests);
  int Recv(void *buf, int count, MPI_Datatype datatype, int source,
           std::vector<MPI_Request> &requests);

  // Helping functions to provide current rank and partner tags (MyRankId as
  // member variable to avoid multiple calls of Mpi library)
  int TagForRank(unsigned int const partner);
  void ResetTagsForPartner();
  inline int MyRankId() { return my_rank_id_; }
};

#endif /* COMMUNICATION_MANAGER_H */
