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
#ifndef LEVELSET_INITIALIZER_DEFINITIONS_H
#define LEVELSET_INITIALIZER_DEFINITIONS_H

#include <string>
#include <stdexcept>

#include "utilities/string_operations.h"

/**
 * @brief The LevelsetInitializerType enum gives all levelset initializer types a unique identifier to allow automatic selection of correct levelset
 *        functions. No specific fixed index is required for this enum class, since no indexing is done.
 */
enum class LevelsetInitializerType {
   Functional,
   STL,
   Parametric
};

/**
 * @brief Converts a string to its corresponding levelset initializer type.
 * @param levelset_initializer_type The levelset initializer type.
 * @return Type of the levelset initializer.
 */
inline LevelsetInitializerType StringToLevelsetInitializerType( std::string const& levelset_initializer_type ) {
   // transform string to upper case without spaces
   std::string const levelset_initializer_upper_case( StringOperations::ToUpperCaseWithoutSpaces( levelset_initializer_type ) );
   // switch statements cannot be used with strings
   if( levelset_initializer_upper_case == "FUNCTIONAL" ) {
      return LevelsetInitializerType::Functional;
   } else if( levelset_initializer_upper_case == "STL" ) {
      return LevelsetInitializerType::STL;
   } else if( levelset_initializer_upper_case == "PARAMETRIC" ) {
      return LevelsetInitializerType::Parametric;
   } else {
      throw std::logic_error( "Levelset initializer type " + levelset_initializer_type + " is not known!" );
   }
}

#endif// LEVELSET_INITIALIZER_DEFINITIONS_H