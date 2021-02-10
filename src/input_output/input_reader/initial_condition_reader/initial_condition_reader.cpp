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
#include "input_output/input_reader/initial_condition_reader/initial_condition_reader.h"

#include "utilities/string_operations.h"

/**
 * @brief Reads the initial condition expression string for a given material index.
 * @param material_index Index of the material for which the data is read.
 * @return Expression string for the material initial condition.
 */
std::string InitialConditionReader::ReadMaterialInitialConditions( unsigned int const material_index ) const {
   return DoReadMaterialInitialConditions( material_index );
}

/**
 * @brief Gives the levelset initialization type name for a given levelset index.
 * @param levelset_index The index of the levelset field (index start: 1).
 * @return Identifier of the levelset initialization type name.
 */
LevelsetInitializerType InitialConditionReader::ReadLevelsetInitializerType( unsigned int const levelset_index, LevelsetInitializerType const default_type ) const {
   std::string const initializer_type( DoReadLevelsetInitializerType( levelset_index ) );
   return initializer_type.empty() ? default_type : StringToLevelsetInitializerType( initializer_type );
}

/**
 * @brief Gives the levelset initialization input for a given levelset index.
 * @param levelset_index The index of the levelset field (index start: 1).
 * @return The levelset initializer input string.
 */
std::string InitialConditionReader::ReadLevelsetInitializerInput( unsigned int const levelset_index ) const {
   std::string const initializer_input( DoReadLevelsetInitializerInput( levelset_index ) );
   if( initializer_input.empty() ) {
      throw std::invalid_argument( "The input for the levelset initial condition must be non-empty!" );
   }
   return StringOperations::RemoveSpaces( initializer_input );
}

/**
 * @brief Gives the variables for the parametric initializer.
 * @param levelset_index The index of the levelset (index start: 1).
 * @return The parametric variables.
 */
std::vector<ParametricVariable> InitialConditionReader::ReadParametricLevelsetInitializerVariables( unsigned int const levelset_index ) const {
   std::vector<std::tuple<std::string, double, double, std::uint64_t>> const parametric_variable_data( DoReadParametricLevelsetInitializerVariables( levelset_index ) );
   std::vector<ParametricVariable> parametric_variables( parametric_variable_data.size() );
   std::transform( std::cbegin( parametric_variable_data ), std::cend( parametric_variable_data ), std::begin( parametric_variables ),
                   []( auto const& data_map ) { return ParametricVariable( std::get<0>( data_map ), std::get<1>( data_map ), std::get<2>( data_map ), std::get<3>( data_map ) ); } );
   return parametric_variables;
}

/**
 * @brief Gives the levelset bounding boxes for a given single levelset index.
 * @param levelset_index The index of the levelset (index start: 1).
 * @return Bounding boxes
 */
std::array<double, 3> InitialConditionReader::ReadParametricLevelsetInitializerReferencePoint( unsigned int const levelset_index ) const {
   return DoReadParametricLevelsetInitializerReferencePoint( levelset_index );
}

/**
 * @brief Gives the levelset bounding boxes for a given single levelset index.
 * @param levelset_index The index of the levelset (index start: 1).
 * @return Bounding boxes
 */
std::vector<std::array<double, 6>> InitialConditionReader::ReadLevelsetInitializerBoundingBoxes( unsigned int const levelset_index ) const {
   return DoReadLevelsetInitializerBoundingBoxes( levelset_index );
}
