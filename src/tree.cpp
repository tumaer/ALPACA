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

#include "tree.h"

#include <utility>
#include <algorithm>
#include <mpi.h>

#include "node.h"
#include "boundary_condition/boundary_factory.h"

/**
 * @brief Default constructor. Creates all nodes on level zero which should reside on the mpi rank the tree is on itself.
 * @param setup Reference used to incorporate user settings, in particular the geometry of the domain and the fluids present at initial condition.
 * @param topology Reference to obtain information about the global state of the simulation.
 */
Tree::Tree(const SimulationSetup& setup, const TopologyManager& topology) : setup_(setup), topology_(topology) {

  std::array<double,CC::ICX()> x_coordinates_of_cells_centers;
  std::array<double,CC::ICY()> y_coordinates_of_cells_centers;
  std::array<double,CC::ICZ()> z_coordinates_of_cells_centers;
  std::array<double,3> origin_of_domain_in_node;
  double cell_size;

  for(const auto& iterator : topology_.LocalLevelZeroIds() ) {
    origin_of_domain_in_node = Node::DomainCoordinatesOfId(iterator,setup.LevelZeroBlockSize());
    cell_size = Node::DomainSizeOfId(iterator,setup.LevelZeroBlockSize()) / CC::ICX(); //ICX is always present
    for(unsigned int i = 0; i < CC::ICX(); ++i) {
        // Corrdinate + Shift to Cell Center + Shift for each Cell
        x_coordinates_of_cells_centers[i] = origin_of_domain_in_node[0] + 0.5*cell_size + i*cell_size;
    }
    for(unsigned int i = 0; i < CC::ICY(); ++i) {
        // Corrdinate + Shift to Cell Center + Shift for each Cell
        y_coordinates_of_cells_centers[i] = origin_of_domain_in_node[1] + 0.5*cell_size + i*cell_size;
    }
    for(unsigned int i = 0; i < CC::ICZ(); ++i) {
        // Corrdinate + Shift to Cell Center + Shift for each Cell
        z_coordinates_of_cells_centers[i] = origin_of_domain_in_node[2] + 0.5*cell_size + i*cell_size;
    }
    InsertNode(iterator, setup_.GetInitialFluids(x_coordinates_of_cells_centers,y_coordinates_of_cells_centers,z_coordinates_of_cells_centers));
  }

  nodes_[0].shrink_to_fit();
}

/**
 * @brief Inserts a Node with specified id and materials into the tree structure.
 * @param id The unique id of the node to be created.
 * @param materials Identifiers of the materials present in the node to be created.
 */
void Tree::InsertNode(const std::uint64_t id,const std::vector<MaterialName> materials) {

  std::array<std::shared_ptr<BoundaryCondition>,6> boundaries;

  boundaries[0] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eEast>(topology_,*this,setup_,id);
  boundaries[1] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eWest>(topology_,*this,setup_,id);
  boundaries[2] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eNorth>(topology_,*this,setup_,id);
  boundaries[3] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eSouth>(topology_,*this,setup_,id);
  boundaries[4] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eTop>(topology_,*this,setup_,id);
  boundaries[5] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eBottom>(topology_,*this,setup_,id);

  nodes_[Node::LevelOfNode(id)].push_back(std::make_shared<Node>(id, boundaries, setup_.LevelZeroBlockSize(), materials));
}

/**
 * @brief Returns pointer to the node having the requested id. Throws exception if node is not in this tree,
 *        i.e. must only be called on correct mpi rank.
 * @param id Id to uniquely identify the Node.
 * @return std::shared_ptr<Node> of requested Node if existent. Throws exception otherwise.
 */
std::shared_ptr<Node> Tree::GetNodeWithId(const std::uint64_t id) const {
  unsigned short level = Node::LevelOfNode(id);

  auto result = std::find_if( nodes_[level].begin(), nodes_[level].end(),[&id](const std::shared_ptr<Node> node) {return (node->GetId() == id);});
  if(result == nodes_[level].end())
    throw std::logic_error("Node with requested Id does not exist: "+std::to_string(id));
  else
    return *result;
}

/**
 * @brief Returns a list of all leaf nodes on this rank. $List is in arbitrary order$.
 * @return List of pointers to the leaves in this tree instance.
 */
std::vector<std::shared_ptr<Node>> Tree::Leaves() const {

  std::vector<std::shared_ptr<Node>> leaves;

  for(const auto& node_id_iterator : topology_.GetLeafIds()) {
    //The node is a, leaf we insert if it is on our rank
    if( topology_.NodeIsOnMyRank(node_id_iterator)) {
        leaves.push_back(GetNodeWithId(node_id_iterator)); //We add this leaf
    }
  }

  return leaves;
}

/**
 * @brief Refines given Node in a three dimensional simulation, inserts its eight children into this tree instance and returns
 *        the ids of the created children. $NOT SAFE: Wrong input gives corrupted data/tree structure. Must only be called on leaves$.
 * @param id Id of the leaf which is to be refined, i.e. becomes a parent node.
 * @return List of childrens' ids.
 */
std::array<std::uint64_t, CC::NOC()> Tree::RefineNode(const std::uint64_t id) {
  // Assign id and cut off the shadow levels
  std::vector<std::uint64_t> ids_children;
  unsigned short level = Node::LevelOfNode(id);

  if(level < nodes_.size() - 1) {
    ids_children = Node::IdsOfChildren(id);

    //insert only the first NOC children, i.e. 2 for 1D, 4 for 2D, and 8 for 3D
    for (unsigned int i = 0; i < CC::NOC(); ++i) {
        InsertNode(ids_children[i], GetNodeWithId(id)->GetFluids());
    }
  } else {
    throw std::invalid_argument("Nodes on this level cannot be refeined further");
  }

  // Saving Memory
  for(auto& level : nodes_) {
    level.shrink_to_fit();
  }

  //select only the first CC::NOC() children to return, i.e. 2 for 1D, 4 for 2D, and 8 for 3D
  std::array<std::uint64_t, CC::NOC()> ids_return;
  std::copy_n(ids_children.begin(), CC::NOC(), ids_return.begin());

  return ids_return;
}

/**
 * @brief Removes the node with the given id. $Not Safe: Throws exception if requested node does not exist in this instance$.
 * @param id The identifier of the node to be removed.
 */
void Tree::RemoveNodeWithId(std::uint64_t const id) {

  unsigned short level = Node::LevelOfNode(id);
  auto position = std::find_if( nodes_[level].begin(), nodes_[level].end(),[&id](const std::shared_ptr<Node> node) {return (node->GetId() == id);});

  if(position == nodes_[level].end()) {
    throw std::invalid_argument("Node to remove does not exist");
  } else {
    nodes_[level].erase(position);
  }

}

/**
 * @brief Allows creation of nodes in this instance from the outside. Should therefore be used with extreme caution.
 *        Needed e.g. for load balancing operations.
 * @param id Unique identifier of the node to be created.
 * @param block The fluid data to be contained in the new node.
 */
void Tree::CreateNode(const std::uint64_t id, Block& block) {

  std::array<std::shared_ptr<BoundaryCondition>,6> bcs;

  bcs[0] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eEast>(topology_,*this,setup_,id);
  bcs[1] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eWest>(topology_,*this,setup_,id);
  bcs[2] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eNorth>(topology_,*this,setup_,id);
  bcs[3] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eSouth>(topology_,*this,setup_,id);
  bcs[4] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eTop>(topology_,*this,setup_,id);
  bcs[5] = BoundaryFactory::DetemineBoundaryType<BoundaryLocation::eBottom>(topology_,*this,setup_,id);

  nodes_[Node::LevelOfNode(id)].push_back(std::make_shared<Node>(id, bcs, setup_.LevelZeroBlockSize(),block));
}

/**
 * @brief Gives a list of all nodes in this instance on the specified level. List is in arbitrary order.
 * @param level The level of interest.
 * @return List of nodes.
 */
const std::vector<std::shared_ptr<Node>>& Tree::NodesOnLevel(const unsigned int level) const {
  if(level > nodes_.size()) {
    throw std::invalid_argument("Requested Level does not exist");
  }

  return nodes_[level];
}

/**
 * @brief Gives a list of all leaves on the specified level. List is in arbitrary order.
 * @param level The level of interest.
 * @return List of leaves.
 */
std::vector<std::shared_ptr<Node>> Tree::LeavesOnLevel(const unsigned int level) const {
  const std::vector<std::shared_ptr<Node>>& leaves = Leaves();
  std::vector<std::shared_ptr<Node>> results;
  std::copy_if(leaves.begin(),leaves.end(),std::back_inserter(results),[&level](const std::shared_ptr<Node>& node){return Node::LevelOfNode(node->GetId()) == level;});
  return results;
}
