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
#ifndef INDEX_TRANSFORMATIONS_H
#define INDEX_TRANSFORMATIONS_H

#include "user_specifications/compile_time_constants.h"

/**
 * @brief This namespace provides functionality to transform buffer indices. "BIT" means "Buffer Index Transformation".
 */
namespace BIT {

   /**
 * @brief Transforms an x-index in an buffer of size "total cells" to an index in an buffer of size "internal cells".
 *        "T2I" means "Total Cell to Inner Cell".
 * @param i The x-index in the buffer of size "total cells".
 * @return The x-index in the buffer of size "inner cells".
 */
   constexpr unsigned int T2IX( unsigned int i ) {
      return i - CC::FICX();
   }

   /**
 * @brief Transforms an y-index in an buffer of size "total cells" to an index in an buffer of size "internal cells".
 *        "T2I" means "Total Cell to Inner Cell".
 * @param j The y-index in the buffer of size "total cells".
 * @return The y-index in the buffer of size "inner cells".
 */
   constexpr unsigned int T2IY( unsigned int j ) {
      return CC::DIM() != Dimension::One ? j - CC::FICY() : 0;
   }

   /**
 * @brief Transforms an z-index in an buffer of size "total cells" to an index in an buffer of size "internal cells".
 *        "T2I" means "Total Cell to Inner Cell".
 * @param k The z-index in the buffer of size "total cells".
 * @return The z-index in the buffer of size "inner cells".
 */
   constexpr unsigned int T2IZ( unsigned int k ) {
      return CC::DIM() == Dimension::Three ? k - CC::FICZ() : 0;
   }

   /**
 * @brief Transforms an x-index in an buffer of size "inner cells" to an index in an buffer of size "total cells".
 *        "I2T" means "Inner Cell to Total Cell".
 * @param i The x-index in the buffer of size "inner cells".
 * @return The x-index in the buffer of size "total cells".
 */
   constexpr unsigned int I2TX( unsigned int i ) {
      return i + CC::FICX();
   }

   /**
 * @brief Transforms an y-index in an buffer of size "inner cells" to an index in an buffer of size "total cells".
 *        "I2T" means "Inner Cell to Total Cell".
 * @param j The y-index in the buffer of size "inner cells".
 * @return The y-index in the buffer of size "total cells".
 */
   constexpr unsigned int I2TY( unsigned int j ) {
      return CC::DIM() != Dimension::One ? j + CC::FICY() : 0;
   }

   /**
 * @brief Transforms an z-index in an buffer of size "inner cells" to an index in an buffer of size "total cells".
 *        "I2T" means "Inner Cell to Total Cell".
 * @param k The z-index in the buffer of size "inner cells".
 * @return The z-index in the buffer of size "total cells".
 */
   constexpr unsigned int I2TZ( unsigned int k ) {
      return CC::DIM() == Dimension::Three ? k + CC::FICZ() : 0;
   }

   /**
 * @brief Transforms an x-index in an buffer of size "total cells" to an index in an flux buffer.
 *        "T2F" means "Total Cell to Flux Array".
 * @param i The x-index in the buffer of size "total cells".
 * @return The x-index in the flux buffer.
 */
   constexpr unsigned int T2FX( unsigned int i ) {
      return i - ( CC::FICX() - 1 );
   }

   /**
 * @brief Transforms an y-index in an buffer of size "total cells" to an index in an flux buffer.
 *        "T2F" means "Total Cell to Flux Array".
 * @param j The y-index in the buffer of size "total cells".
 * @return The y-index in the flux buffer.
 */
   constexpr unsigned int T2FY( unsigned int j ) {
      return j - ( ( CC::DIM() != Dimension::One ) ? ( CC::FICY() - 1 ) : ( -1 ) );
   }

   /**
 * @brief Transforms an z-index in an buffer of size "total cells" to an index in an flux buffer.
 *        "T2F" means "Total Cell to Flux Array".
 * @param k The z-index in the buffer of size "total cells".
 * @return The z-index in the flux buffer.
 */
   constexpr unsigned int T2FZ( unsigned int k ) {
      return k - ( ( CC::DIM() == Dimension::Three ) ? ( CC::FICZ() - 1 ) : ( -1 ) );
   }

}// namespace BIT

#endif// INDEX_TRANSFORMATIONS_H
