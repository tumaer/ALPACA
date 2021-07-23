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
#ifndef FIELD_BUFFER_H
#define FIELD_BUFFER_H

#include "user_specifications/compile_time_constants.h"
#include "block_definitions/field_material_definitions.h"
#include "block_definitions/field_interface_definitions.h"

/**
 * @brief Bundles buffers for material fields of a certain type to have them contiguous in memory and allows accessing each field separately.
 * @tparam N Number of fields.
 * @tparam FieldEnum Enumeration type allowing to access the material fields.
 * @tparam int(*const FieldToIndex)(FieldEnum) Function converting the field enumeration to an index in the range [0;N).
 */
template<std::size_t N, typename FieldEnum, unsigned int ( *const FieldToIndex )( FieldEnum )>
struct FieldBuffer {
   std::array<double[CC::TCX()][CC::TCY()][CC::TCZ()], N> Fields;

   /**
    * @brief Access the buffer corresponding to field f.
    */
   auto operator[]( FieldEnum const f ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
      return Fields[FieldToIndex( f )];
   }

   /**
    * @brief Access the buffer corresponding to field f. Const overload.
    */
   auto operator[]( FieldEnum const f ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
      return Fields[FieldToIndex( f )];
   }

   /**
    * @brief Access the buffer at the given index.
    */
   auto operator[]( unsigned short const index ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
      return Fields[index];
   }

   /**
    * @brief Access the buffer at the given index. Const overload.
    */
   auto operator[]( unsigned short const index ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
      return Fields[index];
   }

   /**
    * @brief Allows to get and set the field values in a certain material cell.
    */
   struct CellView {
      FieldBuffer& buffer_;
      unsigned int const i, j, k;

      double& operator[]( FieldEnum const f ) {
         return buffer_[f][i][j][k];
      }

      double& operator[]( unsigned short const index ) {
         return buffer_[index][i][j][k];
      }
   };

   /**
    * @brief Allows to get the field values in a certain material cell. Read-only!
    */
   struct CellViewConst {
      FieldBuffer const& buffer_;
      unsigned int const i, j, k;

      double operator[]( FieldEnum const f ) const {
         return buffer_[f][i][j][k];
      }

      double operator[]( unsigned short const index ) const {
         return buffer_[index][i][j][k];
      }
   };

   /**
    * @brief Get a non-modifiable object representing the material cell at the given indices.
    * @param i, j, k Material cell index.
    * @return A constant view of the field values in the cell.
    */
   auto GetCellView( unsigned int const i, unsigned int const j, unsigned int const k ) const {
      return CellViewConst{ *this, i, j, k };
   }

   /**
    * @brief Get a modifiable object representing the material cell at the given indices.
    * @param i, j, k Material cell index.
    * @return A handle to the field values in the cell.
    */
   auto GetCellView( unsigned int const i, unsigned int const j, unsigned int const k ) {
      return CellView{ *this, i, j, k };
   }

   /**
    * @brief Gives the number of buffers contained in the field.
    * @return Number of fields contained in the buffer.
    */
   static constexpr std::size_t GetNumberOfFields() {
      return N;
   }
};

/**
 * @brief Bundles the conservative values to have them contiguous in memory.
 */
using Conservatives = FieldBuffer<MF::ANOE(), Equation, ETI>;
// Check Memory Layout at compile time for safe MPI sending (Ensures Compiler did not pad the struct)
static_assert( sizeof( Conservatives ) == MF::ANOE() * CC::TCX() * CC::TCY() * CC::TCZ() * sizeof( double ),
               "Conservative Struct is not contiguous in Memory" );

/**
 * @brief Bundles the prime state values to have them contiguous in memory.
 */
using PrimeStates = FieldBuffer<MF::ANOP(), PrimeState, PTI>;
// Check Memory Layout at compile time for safe MPI sending (Ensures Compiler did not pad the using)
static_assert( sizeof( PrimeStates ) == MF::ANOP() * CC::TCX() * CC::TCY() * CC::TCZ() * sizeof( double ),
               "Prime State Struct is not contiguous in Memory" );

/**
 * @brief Bundles the parameters to have them contiguous in memory.
 */
using Parameters = FieldBuffer<MF::ANOPA(), Parameter, PTI>;
static_assert( sizeof( Parameters ) == MF::ANOPA() * CC::TCX() * CC::TCY() * CC::TCZ() * sizeof( double ),
               "Parameters Struct is not contiguous in Memory" );

/**
 * @brief Bundles the interface descriptions to have them contiguous in memory.
 */
using InterfaceDescriptions = FieldBuffer<IF::ANOD(), InterfaceDescription, IDTI>;
static_assert( sizeof( InterfaceDescriptions ) == IF::ANOD() * CC::TCX() * CC::TCY() * CC::TCZ() * sizeof( double ),
               "InterfaceDescription Struct is not contiguous in Memory" );

/**
 * @brief Bundles the interface states to have them contiguous in memory.
 */
using InterfaceStates = FieldBuffer<IF::ANOS(), InterfaceState, ISTI>;
static_assert( sizeof( InterfaceStates ) == IF::ANOS() * CC::TCX() * CC::TCY() * CC::TCZ() * sizeof( double ),
               "InterfaceState Struct is not contiguous in Memory" );

/**
 * @brief Bundles the interface parameters to have them contiguous in memory.
 */
using InterfaceParameters = FieldBuffer<IF::ANOPA(), InterfaceParameter, IPTI>;
static_assert( sizeof( InterfaceParameters ) == IF::ANOPA() * CC::TCX() * CC::TCY() * CC::TCZ() * sizeof( double ),
               "InterfaceParameter Struct is not contiguous in Memory" );

#endif// FIELD_BUFFER_H
