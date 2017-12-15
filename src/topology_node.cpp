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

#include "topology_node.h"
#include "node.h"

/**
 * @brief Standard constructor for a TopologyNode. Initializes the node as a leaf without children and without assigning a future rank.
 * @param id The id of the node to be created.
 * @param rank Rank holding the node to be created.
 */
TopologyNode::TopologyNode(const std::uint64_t id, const int rank) :
    mpi_tag_ub_(ForwardTagUb()),
    unique_id_(id),
    current_rank_(rank),
    future_rank_(-1),
    children_{},
    is_leaf_(true),
    weight_(1), //Value one is fine - real value needs to be knwon only before load balancing
    tag_(-1) {} //Default "-1" ensures tag is reset, otherwiese MPI_error is trifggerd.

/**
 * @brief Gives the MPI_TAG_UB. Avoids handle creation, e.g. for const members in initializer list.
 * @return MPI_TAG_UB
 */
int TopologyNode::ForwardTagUb() const {
  int *tag_ub;
  int flag;
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);
  return *tag_ub;
}

/**
 * @brief Traces the child with the given id and triggers refinement.
 * @param id The id of the node that is to be refined.
 */
void TopologyNode::Refine(const std::uint64_t id) {
  if (unique_id_ == id) {
      if (is_leaf_) {
          Refine();
      }
  } else {
    GetChildWithId(id).Refine(id);
  }
}

/**
 * @brief Refines the current leaf, i.e. inserts children as leaves and changes status away from leaf.
 */
void TopologyNode::Refine() {
  is_leaf_ = false;
  std::uint64_t working_id = unique_id_ << 3;
  for (int i = 0; i < CC::NOC(); i++) {
      children_.emplace_back(TopologyNode(working_id + i, current_rank_));
  }
}

/**
 * @brief Coarses the node with the given id. I.e. removes all its descendants and sets its status to leaf.
 * @param id The id of the node to be coarsened.
 */
void TopologyNode::Coarse(const std::uint64_t id) {
  if (unique_id_ == id) {
    for(auto& child : children_) {
        child.Coarse(child.Id());
    }
    is_leaf_ = true;
    children_.clear();
  } else {
    if (!is_leaf_) { //coarsening might happen on more than one level, do nothing if already coarsened
        GetChildWithId(id).Coarse(id);
    }
  }
}

/**
 * @brief Gives the node count, i.e. how many nodes are descending from this one.
 * @return If no children are present "1", otherwise 1 + descendants count.
 */
int TopologyNode::Count() const {
  int count = 1;
  for (auto& child : children_) {
      count += child.Count();
  }
  return count;
}

/**
 * @brief Gives the depth of the children, i.e. the number of generations in the tree aka. the number of refinement levels.
 * @return Depth.
 * @note The tree depth starts with 1 by definition.
 */
unsigned int TopologyNode::GetDepth() const {
  if(is_leaf_) {
    return 1;
  } else {
    std::vector<unsigned int> child_depth;
    child_depth.reserve(CC::NOC());
    for(const TopologyNode& child : children_) {
        child_depth.emplace_back(child.GetDepth());
    }
    return 1 + *std::max_element(child_depth.begin(),child_depth.end());
  }
}

/**
 * @brief Gives the tag used for the node with the given id.
 * @param id The id of the node of interest.
 * @return The tag.
 */
int TopologyNode::TagForId(const std::uint64_t id) const {
  if(unique_id_ == id) {
    return tag_;
  } else {
    return GetChildWithId(id).TagForId(id);
  }
}

/**
 * @brief Collects the ids of all nodes on the given level.
 * @param level The level of interest
 * @param ids Indirect return parameter
 */
void TopologyNode::IdsOnLevel(const unsigned int level,std::vector<std::uint64_t>& ids) const {
  if(Node::LevelOfNode(unique_id_) == level) {
    ids.emplace_back(unique_id_);
  } else {
    for(const TopologyNode& child : children_) {
        child.IdsOnLevel(level,ids);
    }
  }
}

/**
 * @brief Gives the rank of the node with the given id.
 * @param id The id of the node of interest.
 * @return Current rank of the node.
 */
int TopologyNode::GetRank(const std::uint64_t id) const {
  if (unique_id_ == id)
      return current_rank_;
  else
      return GetChildWithId(id).GetRank(id);
}

/**
 * @brief Gives a list of all leaves in this TopologyNode (tree).
 * @param leaves Indirect return parameter: The found leaves.
 * @return The number of leaves.
 */
int TopologyNode::GetLeafIds(std::vector<std::uint64_t>& leaves) const {
  if (is_leaf_) {
      leaves.push_back(unique_id_);
      return 1;
  } else {
      int sum = 0;
      for (TopologyNode child : children_) {
          sum += child.GetLeafIds(leaves);
      }
      return sum;
  }
}

/**
 * @brief Checks whether or not the node with the given id exists in the tree.
 * @param id The id of interest.
 * @return True if the node exists, false otherwise (If all children respond with false).
 */
bool TopologyNode::NodeExists(const std::uint64_t id) const {
  if (unique_id_ == id)
      return true;
  else {
      if (is_leaf_) {
          return false;
      }
      return GetChildWithId(id).NodeExists(id);
  }
}

/**
 * @brief Gives a reference to the node with the given id.
 * @param id The id of the desired node.
 * @return The node with given id.
 * @note Does not directly return the correct value, gives a value for recursive search. I.e. if a grand-child is searched this function returns the mother of the grandchild.
 */
TopologyNode& TopologyNode::GetChildWithId(const std::uint64_t id) {
  // We have to search through the ancestor of the searched node until we find the chidl of the current one (recursive calls, remmber...)
  std::uint64_t id_on_child_level = id;
  while(Node::LevelOfNode(id_on_child_level) > (Node::LevelOfNode(unique_id_) + 1) ) {
    id_on_child_level = Node::ParentIdOfNode(id_on_child_level);
  }
  return children_[Node::PositionOfNodeAmongSibilings(id_on_child_level)];
}

/**
 * @brief Const variant of GetChildWithId(const std::uint64_t id), see there for details.
 */
const TopologyNode& TopologyNode::GetChildWithId(const std::uint64_t id) const {
  // We have to search through the ancestor of the searched node until we find the chidl of the current one (recursive calls, remmber...)
  std::uint64_t id_on_child_level = id;
  while(Node::LevelOfNode(id_on_child_level) > (Node::LevelOfNode(unique_id_) + 1) ) {
    id_on_child_level = Node::ParentIdOfNode(id_on_child_level);
  }
  return children_[Node::PositionOfNodeAmongSibilings(id_on_child_level)];
}

/**
 * @brief Ensures that the load (sum of weights) is evenly distributed between the ranks
 *        Assigns the balanced topology to the future_rank_ member, which is then to be used in the actual redistribution in the LoadBalancing routine.
 * @return The rank to which this node should be moved.
 * @note Returns the modified member as of usage in recursive call.
 */
int TopologyNode::BalanceTargetRanks() {
  if (is_leaf_) {
    if (future_rank_ >= 0)
        return future_rank_;
    else
        throw std::logic_error("TopologyNode::balance: Rank for leaves must be set!");
  }else {
    std::vector<int> child_ranks;
    child_ranks.reserve(CC::NOC());
    for(TopologyNode& child : children_) {
        child_ranks.push_back(child.BalanceTargetRanks());
    }
    //Creates a count list, i .e child_on_rank_count[pos] = 4 tells that four children are on the rank on which children_[pos] resides on.
    std::vector<unsigned int> child_on_rank_count(CC::NOC(),0);
    unsigned int last_unique_rank_position = 0;

    for (unsigned int i = 0; i < children_.size(); ++i) {
        for (unsigned int j = 0; j <= i; j++) {
            if (child_ranks[i] == child_ranks[j]) {
                child_on_rank_count[j]++;
                if (j > last_unique_rank_position) {
                    last_unique_rank_position = j;
                }
                break;
              }
          }
      }
      // We retrive the rank on which of the highest count.
      future_rank_ = child_ranks[std::distance(child_on_rank_count.begin(),std::max_element(child_on_rank_count.begin(),child_on_rank_count.end()))];
  }
  return future_rank_;
}

/**
 * @brief Assigns the future ranks according to a homogeneous weight distribution to the leaf TopologyNodes.
 * @param count_rank_map Information bundle mapping the current load "count" to the ranks.
 * @param level The level where the future rank is reassigned
 */
void TopologyNode::SetTargetRankForLeaf(std::vector<std::vector<std::tuple<unsigned int, int>>>& count_rank_map, const unsigned int level) {

  // tupel is [level](<count, rank>)
  // Runs through list, if a rank is already full it is taken out of the list and the following nodes are assigned to the next element with enough room and so on.
  if(is_leaf_) {
    bool assigned = false;
    // loop until this node is assigned to a rank (hence, until a rank with enough room for this node is found)
    while (!assigned) {
        if(std::get<0>(count_rank_map[level].back()) >= weight_) {
            // there's enough room for the current node on this rank, so assign it
            std::get<0>(count_rank_map[level].back()) -= weight_;
            future_rank_ = std::get<1>(count_rank_map[level].back());
            assigned = true;
        } else {
            // there's not enough room for the current node, so take the (possible) overhead and carry it over to the next rank
            // (assignment is tried again in next loop)
            int overhead = std::get<0>(count_rank_map[level].back());
            count_rank_map[level].pop_back();
            std::get<0>(count_rank_map[level].back()) += overhead;
        }
    }
  } else {
    for(unsigned int i=0;i<children_.size();i++) {
        children_[i].SetTargetRankForLeaf(count_rank_map, level+1);
    }
  }

}

/**
 * @brief Overloaded Version of SetTargetRankForLeaf(std::vector<std::vector<std::tuple<int, int>>>&, const unsigned int level)
 *        using a Hilbert Curve traversal for the load definition. See details in other function.
 * @param HilbertPosition The position
 */
void TopologyNode::SetTargetRankForLeaf(std::vector<std::vector<std::tuple<unsigned int, int>>>& count_rank_map, HilbertPosition position, int level){

  // tupel is [level](<count, rank>)
  // Runs through list, if a rank is already full it is taken out of the list and the following nodes are assigned to the next element with enough room and so on.
  if(is_leaf_) {
    bool assigned = false;
    // loop until this node is assigned to a rank (hence, until a rank with enough room for this node is found)
    while (!assigned) {
        if(std::get<0>(count_rank_map[level].back()) >= weight_) {
            // there's enough room for the current node on this rank, so assign it
            std::get<0>(count_rank_map[level].back()) -= weight_;
            future_rank_ = std::get<1>(count_rank_map[level].back());
            assigned = true;
        } else {
            // there's not enough room for the current node, so take the (possible) overhead and carry it over to the next rank
            // (assignment is tried again in next loop)
            int overhead = std::get<0>(count_rank_map[level].back());
            count_rank_map[level].pop_back();
            std::get<0>(count_rank_map[level].back()) += overhead;
        }
    }
  } else {
    const std::array<HilbertPosition,8> replacement = SpaceFillingCurves::GetHilbertReplacement(position);
    const std::array<unsigned short,8> order = SpaceFillingCurves::GetHilbertOrder(position);
    for(unsigned int i=0; i < children_.size(); i++) {
        children_[order[i]].SetTargetRankForLeaf(count_rank_map, replacement[i], level+1);
    }
  }
}


/**
 * @brief Assigns the given rank to the node with the given id.
 * @param id The id of the node to be assigned to the given rank.
 * @param rank The rank the node is assigned to.
 */
void TopologyNode::SetCurrentRankOfLeaf(const std::uint64_t id, const int rank) {
  if (id == unique_id_) {
      if (is_leaf_)
          current_rank_ = rank;
      else
          throw std::logic_error("Avoid calling SetRankOfLeaf for non leaves!");
  } else {
      GetChildWithId(id).SetCurrentRankOfLeaf(id, rank);
  }
}

/**
 * @brief Gives a list of all nodes which are to be moved from one MPI rank to another in order to achieve a better load balance.
 * @param my_rank The rank executing this function $Currently not used$
 * @param ids_current_future_rank_map Indirect return parameter, giving a list of all nodes that need to be shifted fomr one mpi rank to another.
 *        The entries are as the name suggest id - current rank - future rank.
 * @note My rank input is currently not used, but needed in future (planned) optimisations.
 */
void TopologyNode::ListUnbalancedNodes(const int my_rank, std::vector<std::tuple<std::uint64_t, int, int>>& ids_current_future_rank_map) {
  if (future_rank_ != current_rank_) {
      ids_current_future_rank_map.push_back(std::make_tuple(unique_id_, current_rank_, future_rank_));
      current_rank_ = future_rank_;
  }

  if (!is_leaf_) {
      for (TopologyNode& child : children_) {
          child.ListUnbalancedNodes(my_rank, ids_current_future_rank_map);
      }
  }
}

/**
 * @brief Gives a list holding the weight of the nodes (children of root node) on each level.
 * @param list Indirect return parameter. Lists the weight on each level.
 * @param currentLevel The level on which the current node resides on.
 */
void TopologyNode::ChildWeight(std::vector<int>& list, const unsigned int currentLevel) const {
  if(is_leaf_) {
    list[currentLevel] += weight_;
  } else{
    for(const TopologyNode& node : children_) {
        node.ChildWeight(list, currentLevel+1);
    }
  }
}

/**
 * @brief Assigns the given weight to the node with the given id. The requested node is found via recursive tree search.
 * @param id The id of the node of interest.
 * @param weight The weight to be assigned to the node of interest.
 */
void TopologyNode::SetWeightOfId(const std::uint64_t id, const unsigned int weight){
  if(unique_id_ == id){
    weight_ = weight;
  }else{
    GetChildWithId(id).SetWeightOfId(id, weight);
  }
}

/**
 * @brief Assigns unique tags to the TopologyNode.
 * @param tag The tag counter to be used and changed so other nodes get a different tag.
 */
void TopologyNode::AssignTag(int& tag) {
  tag_ = tag++;
  if(tag > (mpi_tag_ub_ >> 5)) {
    throw std::logic_error("You sir (or madam) have to many nodes - tags cannot be unique");
  }

  if(!is_leaf_) {
    for(TopologyNode& child : children_) {
        child.AssignTag(tag);
    }
  }
}

// TODO Once SpaceFillingCurves is a full class with own .cpp file these definitions have to be moved there
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::Rxyz;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::Rx_y_z;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::Ryzx;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::Ry_z_x;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::Rzxy;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::Rz_x_y;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::R_xy_z;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::R_x_yz;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::R_yz_x;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::R_y_zx;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::R_zx_y;
constexpr std::array<HilbertPosition,8> SpaceFillingCurves::R_z_xy;
constexpr std::array<unsigned short,8> SpaceFillingCurves::Oxyz;
constexpr std::array<unsigned short,8> SpaceFillingCurves::Ox_y_z;
constexpr std::array<unsigned short,8> SpaceFillingCurves::Oyzx;
constexpr std::array<unsigned short,8> SpaceFillingCurves::Oy_z_x;
constexpr std::array<unsigned short,8> SpaceFillingCurves::Ozxy;
constexpr std::array<unsigned short,8> SpaceFillingCurves::Oz_x_y;
constexpr std::array<unsigned short,8> SpaceFillingCurves::O_xy_z;
constexpr std::array<unsigned short,8> SpaceFillingCurves::O_x_yz;
constexpr std::array<unsigned short,8> SpaceFillingCurves::O_yz_x;
constexpr std::array<unsigned short,8> SpaceFillingCurves::O_y_zx;
constexpr std::array<unsigned short,8> SpaceFillingCurves::O_zx_y;
constexpr std::array<unsigned short,8> SpaceFillingCurves::O_z_curve;
