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
#ifndef BUFFER_OPERATIONS_H
#define BUFFER_OPERATIONS_H

#include "block_definitions/field_buffer.h"
#include <algorithm>

namespace BufferOperations {

   /**
 * @brief Sets all values of a single buffer to a certain given value.
 * @tparam T type of the cell values.
 * @param cells The buffer cells that are set to the value.
 * @param value The value that is set in all cells.
 */
   template<typename T>
   inline void SetSingleBuffer( T ( &cells )[CC::TCX()][CC::TCY()][CC::TCZ()], T const value ) {
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               cells[i][j][k] = value;
            }
         }
      }
   }

   /**
 * @brief Copies all values from a single buffer into another.
 * @tparam T type of the cell values.
 * @param cells_source The buffer holding the cell values that are copied.
 * @param cells_target The buffer holding the cells where the values are copied into.
 */
   template<typename T>
   inline void CopySingleBuffer( T const ( &cells_source )[CC::TCX()][CC::TCY()][CC::TCZ()], T ( &cells_target )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               cells_target[i][j][k] = cells_source[i][j][k];
            }
         }
      }
   }

   /**
 * @brief Swaps all values from a single buffer with another.
 * @tparam T type of the cell values.
 * @param cells_first The buffer holding the cell values of the first buffer that is swapped.
 * @param cells_second The buffer holding the cell values of the second buffer that is swapped with the first buffer.
 */
   template<typename T>
   inline void SwapSingleBuffer( T ( &cells_first )[CC::TCX()][CC::TCY()][CC::TCZ()], T ( &cells_second )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {
      std::swap( cells_first, cells_second );
   }

   /**
 * @brief Copies a all fields from a FieldBuffer to another.
 * @tparam BufferType The type of the buffer that are copied.
 * @param source_buffer The field buffer from which the values are copied.
 * @param target_buffer The field buffer where the values are copied into.
 */
   template<typename BufferType>
   inline void CopyFieldBuffer( BufferType const& source_buffer, BufferType& target_buffer ) {
      for( size_t field_index = 0; field_index < BufferType::GetNumberOfFields(); ++field_index ) {
         double const( &cells_source )[CC::TCX()][CC::TCY()][CC::TCZ()] = source_buffer[field_index];
         double( &cells_target )[CC::TCX()][CC::TCY()][CC::TCZ()]       = target_buffer[field_index];
         CopySingleBuffer( cells_source, cells_target );
      }
   }

   /**
 * @brief Swaps a all fields from a FieldBuffer to another.
 * @tparam BufferType The type of the buffers that are swapped.
 * @param first_buffer The first field buffer that is part of the swap.
 * @param second_buffer The second field buffer that is part of the swap.
 */
   template<typename BufferType>
   inline void SwapFieldBuffer( BufferType& first_buffer, BufferType& second_buffer ) {
      for( size_t field_index = 0; field_index < BufferType::GetNumberOfFields(); ++field_index ) {
         SwapSingleBuffer( first_buffer[field_index], second_buffer[field_index] );
      }
   }

   /**
 * @brief Sets a a certain fixed value into all fields of a field buffer.
 * @tparam BufferType The type of the buffer that should be set.
 * @tparam T type of the cell values.
 * @param buffer The field buffer where the fixed value is set.
 * @param value The value that is set in the field buffer.
 */
   template<typename BufferType, typename T>
   inline void SetFieldBuffer( BufferType& buffer, T const value ) {
      for( size_t field_index = 0; field_index < BufferType::GetNumberOfFields(); ++field_index ) {
         T( &cells )
         [CC::TCX()][CC::TCY()][CC::TCZ()] = buffer[field_index];
         SetSingleBuffer( cells, value );
      }
   }

   /**
 * @brief Sets certain field-dependent fixed values into the fields of a field buffer.
 * @tparam BufferType The type of the buffer that should be set.
 * @tparam T type of the cell values
 * @param buffer The field buffer where the fixed values are set.
 * @param values The values that are set in the field buffer.
 *
 * @note The values must be in the correct order to be set appropriately.
 */
   template<typename BufferType, typename T>
   inline void SetFieldBuffer( BufferType& buffer, std::array<T, BufferType::GetNumberOfFields()> const& values ) {
      for( size_t field_index = 0; field_index < BufferType::GetNumberOfFields(); ++field_index ) {
         T( &cells )
         [CC::TCX()][CC::TCY()][CC::TCZ()] = buffer[field_index];
         SetSingleBuffer( cells, values[field_index] );
      }
   }

}// namespace BufferOperations

namespace BO = BufferOperations;

#endif// BUFFER_OPERATIONS
