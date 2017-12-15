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

#ifndef TOPOLOGY_MANAGER_H
#define TOPOLOGY_MANAGER_H

#include <cstdint>
#include <vector>
#include <array>
#include <set>
#include <memory>
#include <mpi.h>

#include "simulation_setup.h"
#include "topology_node.h"
#include "compile_time_constants.h"

//No way around this Forward declaration
// If you find one - tell me, please!
class Node;

/**
 * @brief The TopologyManager class handles all aspects relevant for MPI (distributed Memory) parallelization. I.e. overview of data-to-rank maps, sending datatypes, etc.
 *        TopologyManager must not trigger changes in local data. It may only be informed about changes in the local trees. In such a case
 *        it updates the global information automatically and spreads the information to the threads.
 * @note TopologyManager does not have a link to any tree object as it is not interested in tree details, but only in "node counts".
 */
class TopologyManager {

  const SimulationSetup& setup_;

  const int rank_id_;
  const int number_of_ranks_;
  const int mpi_tag_ub_;

  std::array<std::vector<int>, CC::AMNL()> id_to_rank_map_;

  std::vector<std::uint64_t> local_refine_list_;
  std::vector<std::uint64_t> local_coarse_list_;

  std::tuple<std::vector<std::uint64_t>,std::vector<int>> local_weight_list_; //List holds id and weight (of this id) $USED ONLY IN MULTIPHASE VERSION$

  std::vector<TopologyNode> forest_; // A collection of (root) trees is a forest

  // Datatype for Jump-Boundaries, children not listed, cannot have a jump at the given location.
  MPI_Datatype east_jump_child_1_, east_jump_child_3_, east_jump_child_5_, east_jump_child_7_;
  MPI_Datatype west_jump_child_0_, west_jump_child_2_, west_jump_child_4_, west_jump_child_6_;
  MPI_Datatype north_jump_child_2_, north_jump_child_3_, north_jump_child_6_, north_jump_child_7_;
  MPI_Datatype south_jump_child_0_, south_jump_child_1_, south_jump_child_4_, south_jump_child_5_;
  MPI_Datatype top_jump_child_4_, top_jump_child_5_, top_jump_child_6_, top_jump_child_7_;
  MPI_Datatype bottom_jump_child_0_, bottom_jump_child_1_, bottom_jump_child_2_, bottom_jump_child_3_;

  // Datatypes for No-Jump-Boundaries
  MPI_Datatype east_boundary_array_, east_boundary_slice_, west_boundary_array_, west_boundary_slice_, north_boundary_array_, north_boundary_slice_;
  MPI_Datatype south_boundary_array_, south_boundary_slice_, top_boundary_array_, top_boundary_slice_, bottom_boundary_array_, bottom_boundary_slice_;

  /**
   * @brief Creates and commits MPI datatypes for in order to exchange boundary subarrays and such in one go as one subarray.
   */
  template<Dimension DIM>
  void CreateDataTypes();
  /**
   * @brief Frees the memory of the used MPI datatypes.
   */
  template<Dimension DIM>
  void FreeTypes();

  int ForwardRankId() const;
  int ForwardCommunicatorSize() const;
  int ForwardTagUb() const;

  int PositionOfNodeInZeroTopology(const std::uint64_t id) const;
  void AssignBalancedLoad();
  void ListNodeToBalance(std::vector<std::tuple<std::uint64_t, int, int>>& ids_current_future_rank_map);
  void AssignTargetRankToLeaves();

  std::vector<int> WeightsOnLevels() const;

  public:
    TopologyManager() = delete;
    TopologyManager(const SimulationSetup &conf);
    ~TopologyManager();

    int GetTagForNode(const uint64_t id) const;
    int GetRankOfNode(const std::uint64_t id) const;

    unsigned int GetCurrentMaximumLevel() const;
    void GenerateTags();

    void UpdateIdList();

    bool NodeExists(const std::uint64_t id) const;
    bool FaceIsJump(const std::uint64_t id, BoundaryLocation const location) const;
    bool NodeIsOnMyRank(const std::uint64_t id)const;

    void RefineNodeWithId(const std::uint64_t id);
    void CoarseNodeWithId(const std::uint64_t parent_id);
    void PrintStatistics(LogWriter& logger);

    std::vector<std::uint64_t> GetLeafIds() const;
    std::vector<std::uint64_t> LeafIdsOnLevel(const unsigned int level) const;
    std::vector<std::uint64_t> DescendantIdsOfNode(std::uint64_t const id) const;

    std::vector<std::tuple<std::uint64_t, int, int>> GetLoadBalancedTopology();

    MPI_Datatype GetDatatypeForJumpBoundary(const std::pair<uint64_t, BoundaryLocation> jump_boundary) const;
    MPI_Datatype NoJumpSliceDataType(const BoundaryLocation location) const;    //Slices are the internal cells to be sent on the respecitve side of a block. E.g. west slice includes the western most internal cells including north, south, top and bottom halos.
    MPI_Datatype NoJumpBoundaryDataType(const BoundaryLocation location) const;

    void MarkNodeWeight(const std::uint64_t id, const unsigned int weight);

    std::vector<std::uint64_t> LocalLevelZeroIds() const;

    std::vector<std::uint64_t> GlobalIdsOnLevel(unsigned int const level) const;
    /**
     * @brief Gives the id of the invoking MPI Rank. $convenience function to substitute MPI_Comm_rank$
     * @return Rank id.
     */
    inline int MyRankId() const {return rank_id_;}

    /**
     * @brief Gives the overall number of MPI ranks currently running the simulation.
     *        $convenience function to substitute MPI_Comm_size$
     * @return Number of Ranks.
     */
    inline int NumberOfRanks() const {return number_of_ranks_;}

    /**
     * @brief Gives the total number values need to be transferred per conservative variable (+buffer copies),
     *        which need to be transferred if the content of a whole block is to be sent from one rank to another
     */
    static constexpr inline int FullBlockSendingSize() {return CC::TCX()*CC::TCY()*CC::TCZ();}

     /**
     * @brief Gives the total number values need to be transferred per conservative variable (+buffer copies),
     *        if the content of all internal cells of a whole block is to be sent from one rank to another
     */
    static constexpr inline int DomainBlockSendingSize() {return CC::ICX()*CC::ICY()*CC::ICZ();}

    /**
     * @brief Gives the total number of values to be transferred if the jump flux buffers of all conservative variables are to be sent form
     *        one rank to another at once.
     * @return The size of jump flux buffer.
     */
    static constexpr inline int JumpBufferSendingSize() {return CC::NoEq()*CC::ICY()*CC::ICZ();}
};

#endif // TOPOLOGY_MANAGER_H
