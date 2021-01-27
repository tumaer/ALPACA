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
#include "input_output/input_reader/material_reader/material_reader.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Converts a given number of indices into a proper string used as an input tag of the material.
 * @param material_indices Indices to be converted.
 * @return Material input string (e.g., "material1_2_3" for indices {1,2,3}).
 */
std::string MaterialReader::MaterialInputTag( std::vector<unsigned int> const& material_indices ) const {
   // Make vector local
   std::vector<unsigned int> indices( material_indices );
   // Sort the vector to get deterministic behavior independent on order in input vector
   std::sort( indices.begin(), indices.end() );
   // Create proper string
   std::string index_string( "material" + std::to_string( indices.front() ) );
   // erase first already added element
   indices.erase( indices.begin() );
   // Add all indices subsequently
   for( unsigned int const& index : indices ) {
      index_string += "_" + std::to_string( index );
   }
   return index_string;
}

/**
 * @brief Gives the checked number of materials used for the simulation.
 * @return number of materials.
 */
unsigned int MaterialReader::ReadNumberOfMaterials() const {
   // read and make consistency check
   int const number_of_materials( DoReadNumberOfMaterials() );
   if( number_of_materials < 0 ) {
      throw std::invalid_argument( "Number of materials must NOT be below zero!" );
   }

   return static_cast<unsigned int>( number_of_materials );
}

/**
 * @brief Gives the material type for a given material index.
 * @param material_index The index of the material (index start: 1).
 * @param type_default The default type to be used if no is given.
 * @return Identifier of the material type.
 */
MaterialType MaterialReader::ReadMaterialType( unsigned int const material_index, MaterialType const default_type ) const {
   std::string const material_type( DoReadMaterialType( material_index ) );
   return material_type.empty() ? default_type : StringToMaterialType( material_type );
}

/**
 * @brief Gives the equation of state name for a given material index.
 * @param material_index The index of the material (index start: 1).
 * @return Identifier of the equation of state name.
 */
EquationOfStateName MaterialReader::ReadEquationOfStateName( unsigned int const material_index ) const {
   return StringToEos( DoReadEquationOfStateName( material_index ) );
}

/**
 * @brief Gives the equation of state data for a given single material index.
 * @param material_index The index of the material (index start: 1).
 * @return Map storing all provided parameters (no check on consistency and validity).
 */
std::unordered_map<std::string, double> MaterialReader::ReadEquationOfStateData( unsigned int const material_index ) const {
   return DoReadEquationOfStateData( material_index );
}

/**
 * @brief Gives the fixed value of a property for a given set of indices (more than one indicates material pairing models).
 * @param material_indices The indices of the materials (index start: 1).
 * @param property Identifier which property should be read.
 * @return fixed value of the given property.
 */
double MaterialReader::ReadFixedValue( std::vector<unsigned int> const& material_indices, MaterialProperty const property ) const {
   return DoReadFixedValue( material_indices, property );
}

/**
 * @brief Gives the parameter model name of a specific property for a given set of indices (more than one indicates material pairing models).
 * @param material_indices The indices of the materials (index start: 1).
 * @param property Identifier which property should be read.
 * @return Identifier of the parameter model name.
 */
MaterialPropertyModelName MaterialReader::ReadModelName( std::vector<unsigned int> const& material_indices, MaterialProperty const property ) const {
   return StringToMaterialPropertyModel( property, DoReadModelName( material_indices, property ) );
}

/**
 * @brief Gives the model data of a specific property for a given set of indices (more than one indicates material pairing models).
 * @param material_indices The indices of the materials (index start: 1).
 * @param property Identifier which property should be read.
 * @return Map storing all provided model parameters (no check on consistency and validity).
 */
std::unordered_map<std::string, double> MaterialReader::ReadModelData( std::vector<unsigned int> const& material_indices, MaterialProperty const property ) const {
   return DoReadModelData( material_indices, property );
}