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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#ifndef NODE_H
#define NODE_H

#include <vector>
#include <unordered_map>
#include <memory>

#include "boundary_condition/boundary_specifications.h"
#include "materials/material_definitions.h"
#include "block_definitions/block.h"
#include "block_definitions/interface_block.h"
#include "enums/interface_tag_definition.h"
#include "topology/id_information.h"

/**
 * @brief Nodes are the members in the tree. A node holds a block for every phase it contains; the Block then holds the material data.
 *        Node is a container that gathers information common for all phases at a given position, as e.g. Boundary Condition types. Every node
 *        has a unique index for identification, in particular with respect to MPI.
 */
class Node {

   double const node_size_;
   double const node_x_coordinate_;
   std::unordered_map<MaterialName, Block> phases_;

   //type std::int8_t due to definition of enum InterfaceTag. Needs to be changed in case the enum type changes.
   std::int8_t interface_tags_[CC::TCX()][CC::TCY()][CC::TCZ()];
   std::int8_t integrated_interface_tags_[CC::TCX()][CC::TCY()][CC::TCZ()];

   std::unique_ptr<InterfaceBlock> interface_block_;

public:
   Node() = delete;
   explicit Node( nid_t const id, double const node_size_on_level_zero, std::vector<MaterialName> const materials,
                  std::int8_t const initial_interface_tag = ITTI( IT::OldCutCell ) );
   explicit Node( nid_t const id, double const node_size_on_level_zero, std::vector<MaterialName> const materials,
                  std::int8_t const ( &initial_interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()], std::unique_ptr<InterfaceBlock> interface_block = nullptr );
   Node( Node const& ) = delete;
   Node& operator=( Node const& ) = delete;
   Node( Node&& )                 = delete;
   Node& operator=( Node&& ) = delete;
   ~Node()                   = default;

   // Functions to return geometry/topology data of the node
   double GetBlockCoordinateX() const;// NH TODO can we somehow (smart) get rid of this?
   double GetBlockSize() const;
   double GetCellSize() const;

   // Functions to get material data of the node
   Block& GetSinglePhase();
   Block const& GetSinglePhase() const;
   Block& GetPhaseByMaterial( MaterialName const material );
   Block const& GetPhaseByMaterial( MaterialName const material ) const;
   MaterialName GetSinglePhaseMaterial() const;
   std::vector<MaterialName> GetMaterials() const;
   std::unordered_map<MaterialName, Block>& GetPhases();
   std::unordered_map<MaterialName, Block> const& GetPhases() const;

   void AddPhase( MaterialName const material );
   void RemovePhase( MaterialName const material );
   bool ContainsMaterial( MaterialName const material ) const;

   InterfaceBlock& GetInterfaceBlock();
   InterfaceBlock const& GetInterfaceBlock() const;
   void SetInterfaceBlock( std::unique_ptr<InterfaceBlock> interface_block = nullptr );
   bool HasLevelset() const;

   std::int8_t GetUniformInterfaceTag() const;
   template<InterfaceDescriptionBufferType C>
   auto GetInterfaceTags() -> std::int8_t ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   template<InterfaceDescriptionBufferType C>
   auto GetInterfaceTags() const -> std::int8_t const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetInterfaceTags( InterfaceDescriptionBufferType const type ) -> std::int8_t ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetInterfaceTags( InterfaceDescriptionBufferType const type ) const -> std::int8_t const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
};

#endif// NODE_H
