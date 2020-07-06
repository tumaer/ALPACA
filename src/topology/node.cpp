/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#include "node.h"
#include <utility>
#include "id_information.h"

/**
 * @brief Constructs a node object holding blocks for all specified materials.
 * @param id The unique id of this node. $CALLERS RESPONSIBILITY THAT IT IS INDEED UNIQUE!$.
 * @param node_size_on_level_zero The size (= size of internal cells) of a node on level zero.
 * @param materials The materials which are present in this node.
 * @param initial_interface_tag Uniform initial interface tag of the node.
 */
Node::Node(std::uint64_t const id, double const node_size_on_level_zero, std::vector<MaterialName> const materials, std::int8_t const initial_interface_tag) :
   node_size_(DomainSizeOfId(id,node_size_on_level_zero)),
   node_x_coordinate_(DomainCoordinatesOfId(id,node_size_)[0])
{
   for(MaterialName const& material : materials) {
      phases_.emplace(std::piecewise_construct, std::make_tuple(material), std::make_tuple());
   }

   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            interface_tags_[i][j][k] = initial_interface_tag;
         }
      }
   }
}

/**
 * @brief Constructs a node object based on already existing fluid data.
 * @param id The unique id of this node. $CALLERS RESPONSIBILITY THAT IT IS INDEED UNIQUE!$
 * @param node_size_on_level_zero The size (= size of internal cells) of a node on level zero.
 * @param initial_interface_tags Buffer containing the full field of interface tags for this node.
 * @param interface_block Interface block that is added to the node.
 */
Node::Node(std::uint64_t const id, double const node_size_on_level_zero, std::vector<MaterialName> const materials,
   std::int8_t const (&initial_interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()], std::unique_ptr<InterfaceBlock> interface_block) :
   node_size_(DomainSizeOfId(id,node_size_on_level_zero)),
   node_x_coordinate_(DomainCoordinatesOfId(id,node_size_)[0]),
   interface_block_(std::move(interface_block))
{
   for(MaterialName const& material : materials) {
      phases_.emplace(std::piecewise_construct, std::make_tuple(material), std::make_tuple());
   }

   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            interface_tags_[i][j][k] = initial_interface_tags[i][j][k];
         }
      }
   }
}

/**
 * @brief Gives the X-coordinate of the coordinate of the node.
 * @return The X-coordinate of the first (most west-south-bottom) cell in the DOMAIN, i.e. not counting Halos.
 */
double Node::GetBlockCoordinateX() const {
   return node_x_coordinate_;
}

/**
 * @brief Gives the length of the internal part of the node.
 * @return Length of the domain in the node, i.e. cell length * number of internal cells.
 */
double Node::GetBlockSize() const {
   return node_size_;
}

/**
 * @brief Gives the length of a single cell in the blocks of this node.
 * @return cell length.
 */
double Node::GetCellSize() const {
   return node_size_ / CC::ICX();
}

/**
 * @brief Returns the material data in a single-phase node.
 * @return The material data bundled in a Block object.
 */
Block& Node::GetSinglePhase() {
#ifndef PERFORMANCE
   if(phases_.size() > 1) {
      throw std::logic_error("Multi-Nodes do not have a single block");
   }
#endif
   return std::get<1>(*phases_.begin());
}

/**
 * @brief Const overload.
 */
Block const& Node::GetSinglePhase() const {
#ifndef PERFORMANCE
   if(phases_.size() > 1) {
      throw std::logic_error("Multi-Nodes do not have a single block");
   }
#endif
   return std::get<1>(*phases_.cbegin());
}

/**
 * @brief Returns the material of a single phase node.
 * @return The material.
 */
MaterialName Node::GetSinglePhaseMaterial() const {
#ifndef PERFORMANCE
   if(phases_.size() > 1) {
      throw std::logic_error("Multi-Nodes do not have a single material");
   }
#endif
   return std::get<0>(*phases_.cbegin());
}

/**
 * @brief Gives the material identifiers for all phases present in this node.
 * @return The materials.
 */
std::vector<MaterialName> Node::GetMaterials() const {
   std::vector<MaterialName> materials;
   materials.reserve(phases_.size());
   for(auto const& phase : phases_) {
      materials.push_back(phase.first);
   }
   return materials;
}

/**
 * @brief Gives the data of the phases present in this node.
 * @return Vector of block data.
 */
std::unordered_map<MaterialName, Block>& Node::GetPhases() {
   return phases_;
}

/**
 * @brief Const overload.
 */
std::unordered_map<MaterialName, Block> const& Node::GetPhases() const {
   return phases_;
}

/**
 * @brief Returns the material data of the respective material.
 * @param material Name of the material for which the block should be returned.
 * @return The material data as bundled in a Block object.
 */
Block& Node::GetPhaseByMaterial(MaterialName const material) {
   return phases_.at(material);
}

/**
 * @brief Const overload.
 */
Block const& Node::GetPhaseByMaterial(MaterialName const material) const {
   return phases_.at(material);
}

/**
 * @brief Adds an empty block for the given material to the phases of this node.
 * @param material .
 */
void Node::AddPhase(MaterialName const material) {
   // does not test if the material is already present
   phases_.emplace(std::piecewise_construct, std::make_tuple(material), std::make_tuple());
}

/**
 * @brief Removes the block for the given material from the phases of this node.
 * @param material .
 */
void Node::RemovePhase(MaterialName const material) {
   phases_.erase(material);
}

/**
 * @brief Indicates whether or not this node contains the specified material.
 * @param material The material identifier to be checked for.
 * @return True if the material exists in this node. False otherwise.
 */
bool Node::ContainsMaterial(MaterialName const material) const {
   if(phases_.find(material) == phases_.end()) { //C++20 provides a contains function ...
      return false;
   } else {
      return true;
   }
}

/**
 * @brief Returns the InterfaceBlock of the node if it exists. Errors otherwise.
 * @return InterfaceBlock of the Node.
 */
InterfaceBlock& Node::GetInterfaceBlock() {
#ifndef PERFORMANCE
   if(interface_block_ == nullptr) {
      throw std::logic_error("Do not request a InterfaceBlock on a Node that does not have one");
   } else {
      return *interface_block_;
   }
#else
   return *interface_block_;
#endif
}

/**
 * @brief Const overload of GetInterfaceBlock. See there for details.
 */
InterfaceBlock const& Node::GetInterfaceBlock() const {
#ifndef PERFORMANCE
   if(interface_block_ == nullptr) {
      throw std::logic_error("Do not request a InterfaceBlock on a Node that does not have one");
   } else {
      return *interface_block_;
   }
#else
   return *interface_block_;
#endif
}

/**
 * @brief Sets the levelset block of this node. If nullptr is given (default) the current levelset block is released.
 * @param interface_block The new levelset block for this node or empty.
 */
void Node::SetInterfaceBlock(std::unique_ptr<InterfaceBlock> interface_block) {
   interface_block_ = std::move(interface_block);
}

/**
 * @brief Returns the interface tag that is present in all cells for single phase nodes. $Must only be called on single nodes$.
 * @return The representative interface tag for this node.
 */
std::int8_t Node::GetUniformInterfaceTag() const {
   if(phases_.size() > 1) {
      throw std::logic_error("Multi-Nodes do not have uniform tags");
   } else {
      return interface_tags_[CC::FICX()][CC::FICY()][CC::FICZ()]; //NH: We just return the corner. This seems to be the safest option.
   }
}

/**
 * @brief Gives the Interface Tag Buffer.
 * @return Interface tag buffer.
 */
auto Node::GetInterfaceTags() -> std::int8_t (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return interface_tags_;
}

/**
 * @brief Const overlaod. See in non-const for details.
 */
auto Node::GetInterfaceTags() const -> std::int8_t const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return interface_tags_;
}

/**
 * @brief Indicates whether or not the node has a InterfaceBlock.
 * @return True if the node has a InterfaceBlock, false otherwise.
 */
bool Node::HasLevelset() const {
   return interface_block_ == nullptr ? false : true;
}

