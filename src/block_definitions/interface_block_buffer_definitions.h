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
#ifndef INTERFACE_BLOCK_BUFFER_DEFINITIONS_H
#define INTERFACE_BLOCK_BUFFER_DEFINITIONS_H

#include "block_definitions/field_buffer.h"

/**
 * @brief Identifier of all buffers in an interface block (to allow single buffer access, e.g., in halo update or extension).
 */
enum class InterfaceBlockBufferType {
   // interface description buffers
   LevelsetBase,
   LevelsetRightHandSide,
   LevelsetReinitialized,
   LevelsetIntegrated,
   VolumeFractionBase,
   VolumeFractionRightHandSide,
   VolumeFractionReinitialized,
   VolumeFractionIntegrated,
   // interface state buffers
   InterfaceStateVelocity,
   InterfaceStatePressurePositive,
   InterfaceStatePressureNegative,
   // interface parameter buffers
   InterfaceParameterSurfaceTensionCoefficient
};

/**
 * @brief Mapping function to map an interface description buffer to the given interface block buffer type.
 * @param description_field The Identifier of the description field (Levelset or VolumeFraction).
 * @param buffer_type The corresponding interface description buffer type (Base, RightHandSide, Reinitialized).
 * @return Appropriate InterfaceBlockBufferType.
 */
static constexpr InterfaceBlockBufferType MapInterfaceDescritpionToInterfaceBlockBufferType( InterfaceDescription const description_field, InterfaceDescriptionBufferType const buffer_type ) {
   switch( buffer_type ) {
      case InterfaceDescriptionBufferType::Base: {
         switch( description_field ) {
            case InterfaceDescription::Levelset: {
               return InterfaceBlockBufferType::LevelsetBase;
            }
            default: {
               return InterfaceBlockBufferType::VolumeFractionBase;
            }
         }
      } break;
      case InterfaceDescriptionBufferType::RightHandSide: {
         switch( description_field ) {
            case InterfaceDescription::Levelset: {
               return InterfaceBlockBufferType::LevelsetRightHandSide;
            }
            default: {
               return InterfaceBlockBufferType::VolumeFractionRightHandSide;
            }
         }
      } break;
      case InterfaceDescriptionBufferType::Integrated: {
         switch( description_field ) {
            case InterfaceDescription::Levelset: {
               return InterfaceBlockBufferType::LevelsetIntegrated;
            }
            default: {
               return InterfaceBlockBufferType::VolumeFractionIntegrated;
            }
         }
      } break;
      // Last possibility (reinitialized)
      default: {
         switch( description_field ) {
            case InterfaceDescription::Levelset: {
               return InterfaceBlockBufferType::LevelsetReinitialized;
            }
            default: {
               return InterfaceBlockBufferType::VolumeFractionReinitialized;
            }
         }
      }
   }
}

/**
 * @brief Mapping function to map an interface state to the given interface block buffer type.
 * @param state The Identifier of the interface state field (Velocity, PressurePositive or PressureNegative).
 * @return Appropriate InterfaceBlockBufferType.
 */
static constexpr InterfaceBlockBufferType MapInterfaceStateToInterfaceBlockBufferType( InterfaceState const state ) {
   switch( state ) {
      case InterfaceState::Velocity: {
         return InterfaceBlockBufferType::InterfaceStateVelocity;
      }
      case InterfaceState::PressurePositive: {
         return InterfaceBlockBufferType::InterfaceStatePressurePositive;
      }
      // Last possibility
      default: {
         return InterfaceBlockBufferType::InterfaceStatePressureNegative;
      }
   }
}

/**
 * @brief Mapping function to map an interface parameter to the given interface block buffer type.
 * @param parameter The Identifier of the interface parameter field (SurfaceTensionCoefficient).
 * @return Appropriate InterfaceBlockBufferType.
 */
static constexpr InterfaceBlockBufferType MapInterfaceParameterToInterfaceBlockBufferType( InterfaceParameter const parameter ) {
   switch( parameter ) {
      default: {
         return InterfaceBlockBufferType::InterfaceParameterSurfaceTensionCoefficient;
      }
   }
}

/**
 * @brief Mapping function to map an interface field and index to the given interface block buffer type.
 * @param field_type The Identifier of the interface field (States, Parameters or Description).
 * @param field_index Index of the given field.
 * @param buffer_type The corresponding interface buffer type (Base, RightHandSide, Reinitialized).
 * @return Appropriate InterfaceBlockBufferType.
 */
static constexpr InterfaceBlockBufferType MapInterfaceFieldToInterfaceBlockBufferType( InterfaceFieldType const field_type, unsigned int const field_index, InterfaceDescriptionBufferType const buffer_type = InterfaceDescriptionBufferType::RightHandSide ) {
   switch( field_type ) {
      case InterfaceFieldType::Description: {
         return MapInterfaceDescritpionToInterfaceBlockBufferType( IF::ITID( field_index ), buffer_type );
      }
      case InterfaceFieldType::Parameters: {
         return MapInterfaceParameterToInterfaceBlockBufferType( IF::ITIP( field_index ) );
      }
      default: {
         return MapInterfaceStateToInterfaceBlockBufferType( IF::ITIS( field_index ) );
      }
   }
}

#endif// INTERFACE_BLOCK_BUFFER_DEFINITIONS_H
