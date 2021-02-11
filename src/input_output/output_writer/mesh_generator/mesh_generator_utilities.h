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
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Utility functions for generating the mesh of the current simulation.
 */
namespace MeshGeneratorUtilities {

   /**
    * @brief Gives the total number of internal vertices of one single block.
    * @return Number of vertices.
    */
   constexpr unsigned int NumberOfInternalVerticesPerBlock() {
      return ( CC::ICX() + 1 ) * ( CC::ICY() + 1 ) * ( CC::ICZ() + 1 );
   }

   /**
    * @brief Gives the total number of cells of one single block.
    * @return Number of vertices.
    */
   constexpr unsigned int NumberOfTotalVerticesPerBlock() {
      return ( CC::TCX() + 1 ) * ( CC::TCY() + 1 ) * ( CC::TCZ() + 1 );
   }

   /**
    * @brief Gives the total number of internal cells of one single block.
    * @return Number of cells.
    */
   constexpr unsigned int NumberOfInternalCellsPerBlock() {
      return CC::ICX() * CC::ICY() * CC::ICZ();
   }

   /**
    * @brief Gives the total number of cells of one single block.
    * @return Number of cells.
    */
   constexpr unsigned int NumberOfTotalCellsPerBlock() {
      return CC::TCX() * CC::TCY() * CC::TCZ();
   }

   /**
    * @brief Gives the cell size for a given block size.
    * @param block_size Size of the block.
    * @return cell_size of the block.
    *
    * @note since only cubic blocks are used, only one direction is required to determine the size. ICX is used since this value is always filled.
    */
   constexpr double CellSizeForBlockSize( double const block_size ) {
      return block_size / double( CC::ICX() );
   }
}// namespace MeshGeneratorUtilities
