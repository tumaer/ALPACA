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

#include "node.h"
#include <bitset>
#include <memory>
#include <stdexcept>


constexpr std::array<std::uint64_t, CC::AMNL()> Node::kHeadbits_;

/**
 * @brief Constructs a node object holding blocks for all specified materials
 * @param id The unique id of this node. $CALLERS RESPONSISBILITY THAT IT IS INDEED UNIQUE!$
 * @param block_boundaries The conditions imposed at edges (halos) of the block.
 * @param block_size_on_level_zero The size (= size of internal cells) of a block on level zero.
 * @param materials The materials which are present in this node
 */
Node::Node(std::uint64_t const id, const std::array<std::shared_ptr<BoundaryCondition>, 6> block_boundaries,const double block_size_on_level_zero, const std::vector<MaterialName> materials) :
    unique_id_(id),
    block_(materials[0]), //2016-11-07 NH TODO Adapt to MultiFluid
    boundaries_(block_boundaries)
{
  std::uint64_t divisor = 1 << (Node::LevelOfNode(id));
  block_size_ = block_size_on_level_zero / double(divisor);
}

/**
 * @brief Constructs a node object based on already existing fluid data.
 * @param id The unique id of this node. $CALLERS RESPONSISBILITY THAT IT IS INDEED UNIQUE!$
 * @param block_boundaries The conditions imposed at edges (halos) of the block.
 * @param block_size_on_level_zero The size (= size of internal cells) of a block on level zero.
 * @param block The fluid data of the material present in this block.
 */
Node::Node(const uint64_t id,const std::array<std::shared_ptr<BoundaryCondition>, 6> block_boundaries,const double block_size_on_level_zero, Block& block) :
    unique_id_(id),
    block_(block),
    boundaries_(block_boundaries)
{
  std::uint64_t divisor = 1 << (Node::LevelOfNode(id));
  block_size_ = block_size_on_level_zero / double(divisor);
}

/**
 * @brief Gives the coordinates of the coordinate system origin of the internal cells of the node.
 * @return The coordinates of the first (most west-south-bottom) cell in the DOMAIN, i.e. not counting Halos.
 * @note The coordinate is the west-south-bottom corner point not the coordinate of the cell center.
 */
std::array<double,3> Node::GetBlockDomainCoordinates() const {

  std::uint64_t level_operator = unique_id_;
  unsigned short last_bit = (level_operator & 0x1);
  double size = block_size_;

  double x = 0;
  double y = 0;
  double z = 0;

  while(level_operator != 10) {
    x += last_bit * size;
    level_operator >>= 1; //Shift by one to the right
    last_bit = (level_operator & 0x1);
    y += last_bit * size;
    level_operator >>= 1; //Shift by one to the right
    last_bit = (level_operator & 0x1);
    z += last_bit * size;
    level_operator >>= 1; //Shift by one to the right
    last_bit = (level_operator & 0x1);
    size *= 2;
  }

  return {x,y,z};
}

/**
 * @brief Gives the length of the internal part of the node.
 * @return Length of the domain in the node, i.e. cell length * number of internal cells.
 */
double Node::GetBlockSize() const {
  return block_size_;
}

/**
 * @brief Gives the length of a single cell in the blocks of this node.
 * @return cell length.
 */
double Node::GetCellSize() const {
  return block_size_ / CC::ICX();
}

/**
 * @brief Gives the level of the node.
 * @return Level of node.
 */
unsigned short Node::GetLevel() const {
  return Node::LevelOfNode(unique_id_);
}

/**
 * @brief Gives the unique id of the node
 * @return Id of node.
 */
std::uint64_t Node::GetId() const {
  return unique_id_;
}

/**
 * @brief Returns the fluid data in the node.
 * @return The fluid data bundled in a Block object.
 */
Block& Node::GetBlock() {
    return block_;
}

/**
 * @brief Gives the boundary at the provided edge.
 * @param location Direction of the edge.
 * @return Reference to the object which modifies the halo cells in this nodes materials.
 */
BoundaryCondition& Node::GetBoundary(const BoundaryLocation location) {

  switch (location) {
    case BoundaryLocation::eEast : {
        return *(boundaries_[0]);
    }
    break;
    case BoundaryLocation::eWest : {
        return *(boundaries_[1]);
    }
    break;
    case BoundaryLocation::eNorth : {
        return *(boundaries_[2]);
    }
    break;
    case BoundaryLocation::eSouth : {
        return *(boundaries_[3]);
    }
    break;
    case BoundaryLocation::eTop : {
      return *(boundaries_[4]);
    }
    break;
    case BoundaryLocation::eBottom : {
      return *(boundaries_[5]);
    }
    break;
    default:
        throw std::invalid_argument("Boundary Location not Found in Node::GetBoundary() - Impossible Error");
    break;
  }
}

/**
 * @brief Returns the material identifiers of the fluids present in the block.
 * @return Material identifers for each fluid.
 */
std::vector<MaterialName> Node::GetFluids() const {
    return {block_.GetMaterial()}; //NH 2016-11-07 TODO Adapt this in Multi Fluid Case
}

/**
 * @brief Calculates the id of the eastern neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of eastern neighbor.
 */
std::uint64_t Node::EastNeigborOfNodeWithId(const std::uint64_t id) {

  unsigned short level = Node::LevelOfNode(id);

  std::bitset<1> indicator = 0;
  std::bitset<64> working_id(id - kHeadbits_[level]);

  for(std::uint8_t i= 0; i < 60; i+=3) {
    working_id[i] = ( (~working_id[i] & ~indicator[0]) | (working_id[i] & indicator[0]) );
    indicator[0] = (working_id[i] | indicator[0]);
  }

  return (working_id.to_ullong() + kHeadbits_[level]);
}

/**
 * @brief Calculates the id of the western neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of western neighbor.
 */
std::uint64_t Node::WestNeigborOfNodeWithId(const std::uint64_t id) {

  unsigned short level = Node::LevelOfNode(id);

  std::bitset<1> indicator = 0;
  std::bitset<64> working_id(id - kHeadbits_[level]);

  for(std::uint8_t i= 0; i < 60; i+=3) {
    working_id[i] = ( (~working_id[i] & ~indicator[0]) | (working_id[i] & indicator[0]) );
    indicator[0] = (~working_id[i] | indicator[0]);
  }

  return (working_id.to_ullong() + kHeadbits_[level]);
}

/**
 * @brief Calculates the id of the northern neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of northern neighbor.
 */
std::uint64_t Node::NorthNeigborOfNodeWithId(const std::uint64_t id) {

  unsigned short level = Node::LevelOfNode(id);

  std::bitset<1> indicator = 0;
  std::bitset<64> working_id(id - kHeadbits_[level]);

  for(std::uint8_t i= 1; i < 60; i+=3) {
    working_id[i] = ( (~working_id[i] & ~indicator[0]) | (working_id[i] & indicator[0]) );
    indicator[0] = (working_id[i] | indicator[0]);
  }

  return (working_id.to_ullong() + kHeadbits_[level]);
}

/**
 * @brief Calculates the id of the southern neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of southern neighbor.
 */
std::uint64_t Node::SouthNeigborOfNodeWithId(const std::uint64_t id) {

  unsigned short level = Node::LevelOfNode(id);

  std::bitset<1> indicator = 0;
  std::bitset<64> working_id(id - kHeadbits_[level]);

  for(std::uint8_t i= 1; i < 60; i+=3) {
    working_id[i] = ( (~working_id[i] & ~indicator[0]) | (working_id[i] & indicator[0]) );
    indicator[0] = (~working_id[i] | indicator[0]);
  }

  return (working_id.to_ullong() + kHeadbits_[level]);
}

/**
 * @brief Calculates the id of the top neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of top neighbor.
 */
std::uint64_t Node::TopNeigborOfNodeWithId(const std::uint64_t id) {

  unsigned short level = Node::LevelOfNode(id);

  std::bitset<1> indicator = 0;
  std::bitset<64> working_id(id - kHeadbits_[level]);

  for(std::uint8_t i= 2; i < 60; i+=3) {
    working_id[i] = ( (~working_id[i] & ~indicator[0]) | (working_id[i] & indicator[0]) );
    indicator[0] = (working_id[i] | indicator[0]);
  }

  return (working_id.to_ullong() + kHeadbits_[level]);
}

/**
 * @brief Calculates the id of the bottom neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @return Id of bottom neighbor.
 */
std::uint64_t Node::BottomNeigborOfNodeWithId(const uint64_t id) {

  unsigned short level = Node::LevelOfNode(id);

  std::bitset<1> indicator = 0;
  std::bitset<64> working_id(id - kHeadbits_[level]);

  for(std::uint8_t i= 2; i < 60; i+=3) {
    working_id[i] = ( (~working_id[i] & ~indicator[0]) | (working_id[i] & indicator[0]) );
    indicator[0] = (~working_id[i] | indicator[0]);
  }

  return (working_id.to_ullong() + kHeadbits_[level]);
}

/**
 * @brief Gives the id of a neighbor at the provided direction.
 * @param id The id of the node whose neighbor is to be found.
 * @param location Direction in which the neighbor is located.
 * @return Id of the neighbor.
 */
std::uint64_t Node::GetNeighborId(const std::uint64_t id, const BoundaryLocation location) {

  switch (location) {
    case BoundaryLocation::eEast:
        return Node::EastNeigborOfNodeWithId(id);
    break;
    case BoundaryLocation::eWest:
        return Node::WestNeigborOfNodeWithId(id);
    break;
    case BoundaryLocation::eNorth:
        return Node::NorthNeigborOfNodeWithId(id);
    break;
    case BoundaryLocation::eSouth:
        return Node::SouthNeigborOfNodeWithId(id);
    break;
    case BoundaryLocation::eTop:
        return Node::TopNeigborOfNodeWithId(id);
    break;
    case BoundaryLocation::eBottom:
        return Node::BottomNeigborOfNodeWithId(id);
    break;
    default:
        throw std::invalid_argument("Boundary Location not Found in Node::GetNeighborId() - Impossible Error");
    break;
  }
}

/**
 * @brief Determines wether a block edge is also an external or domain edge.
 * @param location The direction of the edge under question.
 * @param id The id of the node under investigation.
 * @param setup The simulation settings as provided by the user.
 * @return True if the edge is a domain edge, false otherweise, i.e. internal edge.
 */
bool Node::IsExternalBoundary(const BoundaryLocation location, const uint64_t id, const SimulationSetup& setup) {

  std::bitset<64> input(id - kHeadbits_[Node::LevelOfNode(id)]);

  bool decision = false;

  switch (location) {
    case BoundaryLocation::eWest : {
        std::bitset<64> west_mask(0x249249249249249);
        std::bitset<64> result;
        result = input & west_mask;
        if(!result.to_ullong()) {
            decision = true;
        }
    }
    break;
    case BoundaryLocation::eSouth : {
        std::bitset<64> south_mask(0x492492492492492);
        std::bitset<64> result;
        result = input & south_mask;
        if(!result.to_ullong()) {
            decision = true;
        }
    }
    break;
    case BoundaryLocation::eBottom : {
        std::bitset<64> bottom_mask(0x924924924924924);
        std::bitset<64> result;
        result = input & bottom_mask;
        if(!result.to_ullong()) {
            decision = true;
        }
    }
    break;
    case BoundaryLocation::eEast : {
        std::bitset<1> indicator(1);
        std::uint32_t tows_exponent = 1;
        tows_exponent <<= Node::LevelOfNode(id);
        std::bitset<20> east_mask( (setup.GetLevelZeroBlocksX() * tows_exponent) - 1 );
        for(unsigned int i = 0; i < 20; ++i) {
            indicator[0] = ( ~(input[i*3] ^ east_mask[i]) & indicator[0] );
        }
        if(indicator[0]) {
            decision = true;
        }
    }
    break;
    case BoundaryLocation::eNorth : {
        std::bitset<1> indicator(1);
        std::uint32_t tows_exponent = 1;
        tows_exponent <<= Node::LevelOfNode(id);;
        std::bitset<20> north_mask( (setup.GetLevelZeroBlocksY() * tows_exponent) - 1 );
        for(unsigned int i = 0; i < 20; ++i) {
            indicator[0] = ( ~(input[(i*3)+1] ^ north_mask[i]) & indicator[0] );
        }
        if(indicator[0]) {
            decision = true;
        }
    }
    break;
    case BoundaryLocation::eTop : {
        std::bitset<1> indicator(1);
        std::uint32_t tows_exponent = 1;
        tows_exponent <<= Node::LevelOfNode(id);;
        std::bitset<20> top_mask( (setup.GetLevelZeroBlocksZ() * tows_exponent) - 1 );
        for(unsigned int i = 0; i < 20; ++i) {
            indicator[0] = ( ~(input[(i*3)+2] ^ top_mask[i]) & indicator[0] );
        }
        if(indicator[0]) {
            decision = true;
        }
    }
    break;
    default: {
        throw std::invalid_argument("Boundary Type in Node::IsExternal does not exist");
    }
    break;
  }

  return decision;
}

/**
 * @brief Gives the length of the internal domain within the node of the given id.
 * @param id Id of the node whose size is to be evaluated.
 * @param block_size_on_level_zero The size (= size of internal cells) of a block on level zero.
 * @return Length of the internal domain in the node.
 */
double Node::DomainSizeOfId(const std::uint64_t id, const double block_size_on_level_zero) {
  std::uint64_t divisor = 1 << (Node::LevelOfNode(id));
  return block_size_on_level_zero / double(divisor);
}

/**
 * @brief Gives the level that the given id is associated to. $NOT SAFE, WRONG INPUT RESULTS IN QUITE INCORRECT DATA$
 * @param id Id of node to be evaluated.
 * @return Level of the id.
 */
unsigned int Node::LevelOfNode(const std::uint64_t id) {

  std::uint64_t working_id = (id  >> 21); //Cut Shadows
  unsigned int level = 0;

  while(working_id > 10) {
     working_id >>= 3;
        level++;
  }

  return level;
}
/**
 * @brief Gives the inverse direction. Convenience function.
 * @param location The direction that is to be inverted.
 * @return The inverse direction.
 */
BoundaryLocation Node::OppositeDirection(BoundaryLocation location) {

  switch (location) {
    case BoundaryLocation::eEast:
        return BoundaryLocation::eWest;
    break;
    case BoundaryLocation::eWest:
        return BoundaryLocation::eEast;
    break;
    case BoundaryLocation::eNorth:
        return BoundaryLocation::eSouth;
    break;
    case BoundaryLocation::eSouth:
        return BoundaryLocation::eNorth;
    break;
    case BoundaryLocation::eTop:
        return BoundaryLocation::eBottom;
    break;
    case BoundaryLocation::eBottom:
        return BoundaryLocation::eTop;
    break;
    default:
        throw std::invalid_argument("Error in OppositeDirection - Should never happen");
    break;
  }
}

/**
 * @brief Gives the coordinates of the coordinate system origin of a nodes internal cells.
 * @param id The id of a node for which the coordinates are calculated.
 * @param level_zero_block_size The size (= size of internal cells) of a block on level zero.
 * @return The coordinates of the coordinate origin in the block, i.e. coordinates of corner of first internal cells.
 */
std::array<double, 3> Node::DomainCoordinatesOfId(const uint64_t id, const double level_zero_block_size) {
  std::uint64_t level_operator = id;
  unsigned short last_bit = (level_operator & 0x1);
  double size = Node::DomainSizeOfId(id, level_zero_block_size);

  double x = 0;
  double y = 0;
  double z = 0;

  while(level_operator != 10) {
    x += last_bit * size;
    level_operator >>= 1; //Shift by one to the right
    last_bit = (level_operator & 0x1);
    y += last_bit * size;
    level_operator >>= 1; //Shift by one to the right
    last_bit = (level_operator & 0x1);
    z += last_bit * size;
    level_operator >>= 1; //Shift by one to the right
    last_bit = (level_operator & 0x1);
    size *= 2;
  }

  return {x,y,z};
}
