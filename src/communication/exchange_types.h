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
#ifndef EXCHANGE_TYPES_H
#define EXCHANGE_TYPES_H

#include "user_specifications/compile_time_constants.h"

/**
 * @brief Enum class to give the index of the appropriate Exchanges Types.
 * @note Do not change the underlying type. Used for index mapping.
 */
enum class ExchangeType : unsigned short { Plane = 0,
                                           Stick = 1,
                                           Cube  = 2 };

/**
 * @brief Converts an ExchangeType identifier to a (C++11 standard compliant, i. e. positive) array index. "ETTI = Exchange Type To Index".
 * @param et The exchange type identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<ExchangeType>::type ETTI( ExchangeType const et ) { return static_cast<typename std::underlying_type<ExchangeType>::type>( et ); }

/**
 * @brief Exchange type for a single plane (east, west, north, south, top, bottom).
 */
struct ExchangePlane {
   double plane_[CC::HSSX()][CC::ICY()][CC::ICZ()];
};
static_assert( sizeof( ExchangePlane ) == CC::HSSX() * CC::ICY() * CC::ICZ() * sizeof( double ), "ExchangePlane Struct is not contiguous in Memory" );

/**
 * @brief Exchange type for a single stick representing the intersection between two planes (e.g., north-east).
 */
struct ExchangeStick {
   double stick_[CC::HSSX()][CC::HSSY()][CC::ICZ()];
};
static_assert( sizeof( ExchangeStick ) == CC::HSSX() * CC::HSSY() * CC::ICZ() * sizeof( double ), "ExchangeStick Struct is not contiguous in Memory" );

/**
 * @brief Exchange type for a single cube representing the intersection between three planes (e.g., north-east-bottom).
 */
struct ExchangeCube {
   double cube_[CC::HSSX()][CC::HSSY()][CC::HSSZ()];
};
static_assert( sizeof( ExchangeCube ) == CC::HSSX() * CC::HSSY() * CC::HSSZ() * sizeof( double ), "ExchangeCube Struct is not contiguous in Memory" );

/**
 * @brief Container of the interface tag array. In order to store it in std::vectors. Used for more efficient MPI communication.
 */
struct InterfaceTagBundle {
   std::int8_t interface_tags_[CC::TCX()][CC::TCY()][CC::TCZ()];
};
// Check Memory Layout at compile time for safe MPI sending (Ensures Compiler did not pad the struct)
static_assert( sizeof( InterfaceTagBundle ) == CC::TCX() * CC::TCY() * CC::TCZ() * sizeof( std::int8_t ), "InterfaceTagBundle is not contiguous in Memory" );

#endif// EXCHANGE_TYPES_H
