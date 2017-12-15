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

#ifndef NODE_H
#define NODE_H

#include <cstdint> //64bit ensured ints
#include <array>
#include <vector>
#include <memory> //shared_ptr
#include "simulation_setup.h"
#include "boundary_condition/boundary_specifications.h"
#include "boundary_condition/boundary_condition.h"

#include "block.h"

/**
 * @brief Nodes are the members in the tree. A node holds a block for every phase it contains; the Block then holds the fluid data.
 *        Node is a container that gathers information common for all phases at a given position, as e.g. Boundary Condition types. Every node
 *        has a unique index for identification, in particular with respect to MPI.
 */
class Node {

  // NH: DO NOT TOUCH UNLESS YOU HAVE MY WRITTEN PERMISSION!
  static constexpr std::array<std::uint64_t, CC::AMNL()> kHeadbits_ =
                    {{0x1400000,           0xA000000,          0x50000000,         0x280000000,        0x1400000000,       0xA000000000,
                      0x50000000000,       0x280000000000,     0x1400000000000,    0xA000000000000,    0x50000000000000,   0x280000000000000,
                      0x1400000000000000,  0xA000000000000000 }};

  std::uint64_t const unique_id_;
  double block_size_;
  Block block_;

  const std::array<std::shared_ptr<BoundaryCondition>, 6> boundaries_;

  public:
    Node() = delete;
    Node(std::uint64_t const id, const std::array<std::shared_ptr<BoundaryCondition>, 6> block_boundaries, const double block_size_on_level_zero, const std::vector<MaterialName> materials);
    Node(std::uint64_t const id, const std::array<std::shared_ptr<BoundaryCondition>, 6> block_boundaries, const double block_size_on_level_zero, Block& block);

    std::array<double,3> GetBlockDomainCoordinates() const;

    double GetBlockSize() const;
    double GetCellSize() const;

    unsigned short GetLevel() const;
    std::uint64_t GetId() const;

    Block& GetBlock();

    BoundaryCondition& GetBoundary(const BoundaryLocation location);

    std::vector<MaterialName> GetFluids() const;

    static std::uint64_t EastNeigborOfNodeWithId(const std::uint64_t id);
    static std::uint64_t WestNeigborOfNodeWithId(const std::uint64_t id);
    static std::uint64_t NorthNeigborOfNodeWithId(const uint64_t id);
    static std::uint64_t SouthNeigborOfNodeWithId(const std::uint64_t id);
    static std::uint64_t TopNeigborOfNodeWithId(const std::uint64_t id);
    static std::uint64_t BottomNeigborOfNodeWithId(const std::uint64_t id);

    static std::uint64_t GetNeighborId(const std::uint64_t id, const BoundaryLocation location);

    static bool IsExternalBoundary(BoundaryLocation const location, std::uint64_t const id, const SimulationSetup& setup);

    static double DomainSizeOfId(const std::uint64_t id, const double block_size_on_level_zero);

    static unsigned int LevelOfNode(const uint64_t id);

    static BoundaryLocation OppositeDirection(BoundaryLocation location);

    static std::array<double,3> DomainCoordinatesOfId(const std::uint64_t id, const double level_zero_block_size);
    /**
     * @brief Gives the unique identifer of the parent node. $NO INPUT CHECKS, CALLER MUST ENSURE CORRECT INPUT$
     * @param id Id of the Child
     * @return Id of the Parent
     */
    static inline std::uint64_t ParentIdOfNode(const std::uint64_t id) {return (id >> 3);}

    /**
     * @brief Returns the starting sequence for ids on a level. $NO INPUT CHECKS, CALLER MUST ENSURE CORRECT INPUT$
     * @param level The level of interest. $LEVEL MUST NOT BE OUT OF BOUNDS$
     * @return First Index on the requested level.
     */
    static constexpr inline std::uint64_t HeadbitOnLevel(const unsigned int level) {return kHeadbits_[level];}

    /**
     * @brief Gives the Ids of all eight children of a Node.
     * @param id The Id of the parent node
     * @return Ids of the children in increasing order, i.e. bottom-south-west, bottom-south-east, bottom-north-west, ... ,top-north-east.
     */
     static inline std::vector<std::uint64_t> IdsOfChildren(const std::uint64_t parent_id) {return { {(parent_id << 3),(parent_id << 3)+1
                                                                                            #if DIMENSION > 1
                                                                                                     ,(parent_id << 3)+2,(parent_id << 3)+3
                                                                                            #endif
                                                                                            #if DIMENSION == 3
                                                                                                     ,(parent_id << 3)+4,(parent_id << 3)+5,(parent_id << 3)+6,(parent_id << 3)+7
                                                                                            #endif
                                                                                            }};}


    /**
     * @brief Indicates the position of the Node among its siblings. The position is encoded as integer value with 0 = bottom-south-west to 7 = top-north-east.
     * @param id The id of the node, whose position is to be determined.
     * @return The position among its sibilings
     */
    static inline unsigned short PositionOfNodeAmongSibilings(const std::uint64_t id) {return (id & 0x7);}

};

#endif // NODE_H
