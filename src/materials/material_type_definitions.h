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
#ifndef MATERIAL_TYPE_DEFINITIONS_H
#define MATERIAL_TYPE_DEFINITIONS_H

#include "utilities/string_operations.h"

/**
 * @brief The MaterialType enum gives all types a unique identifier.
 */
enum class MaterialType { Fluid,
                          SolidBoundary
};

/**
 * @brief Converts a string to its corresponding equation of state name.
 * @param eos_name The equation of state name.
 * @return Name of the equation of state.
 */
inline MaterialType StringToMaterialType( std::string const& material_type_name ) {
   // transform string to upper case without spaces
   std::string const type_upper_case( StringOperations::ToUpperCaseWithoutSpaces( material_type_name ) );
   // switch statements cannot be used with strings
   if( type_upper_case == "FLUID" ) {
      return MaterialType::Fluid;
   } else if( type_upper_case == "SOLIDBOUNDARY" ) {
      return MaterialType::SolidBoundary;
   } else {
      throw std::logic_error( "Material type " + type_upper_case + " is not known!" );
   }
}

/**
 * @brief Gives the string for each material type.
 * @param type Material type for which the string should be returned.
 * @return string of material type.
 */
inline std::string MaterialTypeToString( MaterialType const type ) {
   switch( type ) {
      case MaterialType::Fluid: {
         return "Fluid";
      }
      case MaterialType::SolidBoundary: {
         return "Solid boundary";
      }
      default: {
         throw std::logic_error( "Material type not known!" );
      }
   }
}

#endif// MATERIAL_TYPE_DEFINITIONS_H