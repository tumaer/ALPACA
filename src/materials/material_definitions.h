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
#ifndef MATERIAL_DEFINITIONS_H
#define MATERIAL_DEFINITIONS_H

#include <type_traits>
#include <vector>

/**
 * @brief The MaterialName enum gives all materials a unique identifier to allow automatic selection of correct materials as well
 *        as to distinguish treatment in case of contact between materials, e.g. merging or interaction.
 * @note !IMPORTANT! Do NOT change the order of assignment from MaterialOne to index 0, MaterialTwo to index 1, and so fourth. This is used to index the vector
 *       in the material manager class and for assignment of proper material pairings between two material. Adding of additional material by simply increasing numbers
 *       of MaterialOutOfBounds and filling the gap in between.
 */
enum class MaterialName : unsigned short {
   MaterialOne         = 0,
   MaterialTwo         = 1,
   MaterialThree       = 2,
   MaterialFour        = 3,
   MaterialOutOfBounds = 4
};

/**
 * @brief Converts a Material identifier to a (C++11 standard compliant, i. e. positive) sortabel integer. "MTI = Material To Index".
 * @param m The material identifier.
 * @return Integer to be used for sorting.
 */
constexpr std::underlying_type<MaterialName>::type MTI( MaterialName const m ) { return static_cast<typename std::underlying_type<MaterialName>::type>( m ); }

/**
 * @brief Converts a Material index to its corresponding enum name "ITM = Index to Material".
 * @param i The material index.
 * @return Material name for the index.
 */
constexpr MaterialName ITM( typename std::underlying_type<MaterialName>::type i ) { return static_cast<MaterialName>( i ); }

/**
 * @brief The MaterialPairingName enum gives all identifiers for a pairing between two materials.
 * @note !IMPORTANT! Do NOT change the order of assignment from MaterialOneTwo to index 0, MaterialOneThree to index 1, MaterialTwoThree to index 2 and so fourth.
 *       This is used to index the vector in the material manager class and provides a unique mapping function.
 *       If an additionall material is added simply add all pairings in the order MaterialPairing[1:N-1][N] and increase the out of bounds number.
 */
enum class MaterialPairingName : unsigned short {
   MaterialPairingOneTwo      = 0,
   MaterialPairingOneThree    = 1,
   MaterialPairingTwoThree    = 2,
   MaterialPairingOneFour     = 3,
   MaterialPairingTwoFour     = 4,
   MaterialPairingThreeFour   = 5,
   MaterialPairingOutOfBounds = 6
};

/**
 * @brief Converts a MaterialPairing identifier to a (C++11 standard compliant, i. e. positive) sortabel integer. "MPTI = Material Pairing To Index".
 * @param mp The material pairing identifier.
 * @return Integer to be used for sorting.
 */
constexpr std::underlying_type<MaterialPairingName>::type MPTI( MaterialPairingName const mp ) { return static_cast<typename std::underlying_type<MaterialPairingName>::type>( mp ); }

/**
 * @brief Converts a MaterialPairing index to its corresponding enum name "ITMP = Index to Material Pairing".
 * @param i The material pairing index.
 * @return Name of the material pairing.
 */
constexpr MaterialPairingName ITMP( typename std::underlying_type<MaterialPairingName>::type i ) { return static_cast<MaterialPairingName>( i ); }

/**
 * @brief Creates for the given number of materials all pairing indices.
 * @param number_of_materials Number of materials considered.
 * @return vector with pairing indices.
 *
 * @note Do not change the order here, since this reflects the same order as in the enum class of the pairings.
 */
inline std::vector<std::vector<unsigned int>> GetMaterialPairingIndices( unsigned int const number_of_materials ) {
   // Generate all combinations of indices
   std::vector<std::vector<unsigned int>> pairing_indices;
   // Loop structure for creation
   for( unsigned int second_index = 1; second_index < number_of_materials; ++second_index ) {
      for( unsigned int first_index = 0; first_index <= second_index - 1; ++first_index ) {
         pairing_indices.push_back( { first_index + 1, second_index + 1 } );
      }
   }

   return pairing_indices;
}

#endif// MATERIAL_DEFINITIONS_H
