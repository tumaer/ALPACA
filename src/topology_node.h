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

#ifndef TOPOLOGY_NODE_H
#define TOPOLOGY_NODE_H

#include <vector>
#include <tuple>
#include <memory>

#include "space_filling_curves.h"

/**
 * @brief The TopologyNode class organizes the light weight global (over MPI ranks) node information in a tree structure. Allowing the TopologyManager efficient searches.
 */
class TopologyNode {

  const int mpi_tag_ub_;
  const std::uint64_t unique_id_;
  int current_rank_;
  int future_rank_;
  std::vector<TopologyNode> children_;
  bool is_leaf_;
  unsigned int weight_;
  int tag_;

  int ForwardTagUb() const;

  void Refine();

  TopologyNode& GetChildWithId(const std::uint64_t id);
  const TopologyNode& GetChildWithId(const std::uint64_t id) const;

  public:
    TopologyNode() = delete;
    TopologyNode(const std::uint64_t id,const int rank = -1);

    void Refine(const std::uint64_t id);
    void Coarse(const std::uint64_t id);

    int GetRank(const std::uint64_t id) const;
    int Count() const;
    unsigned int GetDepth() const;

    int TagForId(const std::uint64_t id) const;

    void IdsOnLevel(const unsigned int level,std::vector<std::uint64_t>& ids) const;

    int GetLeafIds(std::vector<std::uint64_t>& leaves) const;

    void SetTargetRankForLeaf(std::vector<std::vector<std::tuple<unsigned int,int>>>& count_rank_map, const unsigned int level=0);
    void SetTargetRankForLeaf(std::vector<std::vector<std::tuple<unsigned int,int>>>& count_rank_map, HilbertPosition position, int level=0);

    void ListUnbalancedNodes(const int my_rank, std::vector<std::tuple<std::uint64_t, int, int>>& ids_current_future_rank_map);

    int BalanceTargetRanks();
    void SetCurrentRankOfLeaf(const std::uint64_t id, const int rank);

    bool NodeExists(const std::uint64_t id) const;

    void ChildWeight(std::vector<int>& list, const unsigned int currentLevel=0) const;
    void SetWeightOfId(const std::uint64_t id, const unsigned int weight);

    void AssignTag(int& tag);

    /**
     * @brief Gives the id of this topology node.
     * @return id.
     */
    inline std::uint64_t Id() const {return unique_id_;}

    /**
     * @brief Gives the current rank of the TopologyNode.
     * @return Rank.
     */
    inline int GetRank() const {return current_rank_;}

    inline bool operator==(const std::uint64_t rhs) {return (rhs==unique_id_);}
    inline bool operator==(const TopologyNode& rhs) {return (rhs.Id()== unique_id_);}
};

#endif // TOPOLOGY_NODE_H
