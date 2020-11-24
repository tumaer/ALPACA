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
#ifndef INTERFACE_BLOCK_H
#define INTERFACE_BLOCK_H

#include "user_specifications/compile_time_constants.h"
#include "interface_block_buffer_definitions.h"
#include "block_definitions/field_buffer.h"
#include "block_definitions/field_interface_definitions.h"

/**
 * @brief The InterfaceBlock class holds the interface data, such as the interface describing variables (levelset, volume fraction) and other interfaces
 *        variables that are required for the computations (e.g., interface velocity) Does NOT manipulate the data itself, but provides access to the data.
 */
class InterfaceBlock {

   // buffers for the interface description (different buffer types required for the integration)
   InterfaceDescriptions base_;
   InterfaceDescriptions right_hand_side_;
   InterfaceDescriptions reinitialized_;
   InterfaceDescriptions initial_;
   InterfaceDescriptions integrated_;

   // buffers for the interface states (e.g. interface velocity, negative/positive pressure)
   InterfaceStates states_;

   // buffer to store field-dependent interface parameter (e.g., surface tension coefficient)
   InterfaceParameters parameters_;

public:
   InterfaceBlock() = delete;
   explicit InterfaceBlock( double const levelset_initial );
   explicit InterfaceBlock( double const ( &levelset_initial )[CC::TCX()][CC::TCY()][CC::TCZ()] );
   ~InterfaceBlock()                       = default;
   InterfaceBlock( InterfaceBlock const& ) = delete;
   InterfaceBlock& operator=( InterfaceBlock const& ) = delete;
   InterfaceBlock( InterfaceBlock&& )                 = delete;
   InterfaceBlock& operator=( InterfaceBlock&& ) = delete;

   // Returning general field buffer
   auto GetFieldBuffer( InterfaceFieldType const field_type, unsigned int const field_index, InterfaceDescriptionBufferType const buffer_type = InterfaceDescriptionBufferType::Base ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetFieldBuffer( InterfaceFieldType const field_type, unsigned int const field_index, InterfaceDescriptionBufferType const buffer_type = InterfaceDescriptionBufferType::Base ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];

   // Returning interface description buffers
   auto GetBaseBuffer( InterfaceDescription const interface_description ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetBaseBuffer( InterfaceDescription const interface_description ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetRightHandSideBuffer( InterfaceDescription const interface_description ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetRightHandSideBuffer( InterfaceDescription const interface_description ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetReinitializedBuffer( InterfaceDescription const interface_description ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetReinitializedBuffer( InterfaceDescription const interface_description ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetInitialBuffer( InterfaceDescription const interface_description ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetInitialBuffer( InterfaceDescription const interface_description ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetIntegratedBuffer( InterfaceDescription const interface_description ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetIntegratedBuffer( InterfaceDescription const interface_description ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];

   template<InterfaceDescriptionBufferType C>
   InterfaceDescriptions& GetInterfaceDescriptionBuffer();

   template<InterfaceDescriptionBufferType C>
   InterfaceDescriptions const& GetInterfaceDescriptionBuffer() const;

   InterfaceDescriptions& GetBaseBuffer();
   InterfaceDescriptions const& GetBaseBuffer() const;
   InterfaceDescriptions& GetRightHandSideBuffer();
   InterfaceDescriptions const& GetRightHandSideBuffer() const;
   InterfaceDescriptions& GetReinitializedBuffer();
   InterfaceDescriptions const& GetReinitializedBuffer() const;
   InterfaceDescriptions& GetInitialBuffer();
   InterfaceDescriptions const& GetInitialBuffer() const;
   InterfaceDescriptions& GetIntegratedBuffer();
   InterfaceDescriptions const& GetIntegratedBuffer() const;
   InterfaceDescriptions& GetInterfaceDescriptionBuffer( InterfaceDescriptionBufferType const buffer_type );
   InterfaceDescriptions const& GetInterfaceDescriptionBuffer( InterfaceDescriptionBufferType const buffer_type ) const;

   // Returning state buffers
   auto GetInterfaceStateBuffer( InterfaceState const state_type ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetInterfaceStateBuffer( InterfaceState const state_type ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];

   InterfaceStates& GetInterfaceStateBuffer();
   InterfaceStates const& GetInterfaceStateBuffer() const;

   // returning parameter buffers
   auto GetInterfaceParameterBuffer( InterfaceParameter const parameter_type ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetInterfaceParameterBuffer( InterfaceParameter const parameter_type ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];

   InterfaceParameters& GetInterfaceParameterBuffer();
   InterfaceParameters const& GetInterfaceParameterBuffer() const;

   // returning general interface block buffer
   auto GetBuffer( InterfaceBlockBufferType const buffer_type ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetBuffer( InterfaceBlockBufferType const buffer_type ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()];
};

#endif// INTERFACE_BLOCK_H