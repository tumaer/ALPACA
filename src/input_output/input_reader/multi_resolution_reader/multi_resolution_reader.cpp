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
#include "input_output/input_reader/multi_resolution_reader/multi_resolution_reader.h"

#include <numeric>

#include "user_specifications/compile_time_constants.h"
#include "enums/dimension_definition.h"

/**
 * @brief Gives the checked number of nodes on level zero for a given direction.
 * @param direction Direction identifier for which the number of nodes should be read.
 * @return number of nodes.
 */
unsigned int MultiResolutionReader::ReadNumberOfNodes( Direction const direction ) const {
   // Read the block number and check on consistency
   int const number_of_nodes( DoReadNumberOfNodes( direction ) );
   if( number_of_nodes <= 0 ) {
      throw std::invalid_argument( "At least one block must be present on level zero in each direction!" );
   }
   if( number_of_nodes > 128 ) {
      throw std::invalid_argument( "Block number on level zero must not exceed 128!" );
   }
   return static_cast<unsigned int>( number_of_nodes );
}

/**
 * @brief Gives the checked size of a node on level zero.
 * @return node size on level zero.
 */
double MultiResolutionReader::ReadNodeSizeOnLevelZero() const {
   // Read the node size value and check on consistency
   double const node_size( DoReadNodeSizeOnLevelZero() );
   if( node_size <= 0.0 ) {
      throw std::invalid_argument( "Node size on level zero must be greater zero!" );
   }
   return node_size;
}

/**
 * @brief Gives the checked maximum level used for the simulation.
 * @return maximum level of the simulation.
 */
unsigned int MultiResolutionReader::ReadMaximumLevel() const {
   // Read the value and check on consistency
   int const maximum_level( DoReadMaximumLevel() );
   if( maximum_level > static_cast<int>( CC::AMNL() ) ) {
      throw std::invalid_argument( "Maximum level must NOT be larger than " + std::to_string( CC::AMNL() ) + "!" );
   }
   if( maximum_level < 0 ) {
      throw std::invalid_argument( "Maximum level must NOT be below zero!" );
   }
   return static_cast<unsigned int>( maximum_level );
}

/**
 * @brief Gives the checked reference epsilon value used for the refinement criterion.
 * @return epsilon refernce value.
 */
double MultiResolutionReader::ReadEpsilonReference() const {
   // Read the reference value and check on consistency
   double const reference( DoReadEpsilonReference() );
   if( reference <= 0.0 ) {
      throw std::invalid_argument( "Epsilon reference must be larger than zero!" );
   }
   return reference;
}

/**
 * @brief Gives the checked reference level where the epsilon reference is enforced.
 * @return level of reference.
 */
unsigned int MultiResolutionReader::ReadEpsilonLevelReference() const {
   // Read the reference level and check on consistency
   int const reference_level( DoReadEpsilonLevelReference() );
   if( reference_level < 0 ) {
      throw std::invalid_argument( "Level of epsilon reference must NOT be below zero!" );
   }
   return static_cast<unsigned int>( reference_level );
}