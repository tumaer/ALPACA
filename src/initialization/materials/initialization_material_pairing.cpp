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
#include "initialization/materials/initialization_material_pairing.h"

#include "user_specifications/compile_time_constants.h"
#include "materials/material_property_definitions.h"

namespace Initialization {

   /**
    * @brief Initializes a complete material pairing with the given input reader 
    * @param material_indices Indices of the materials that form the pairing, which should be read from the input file 
    * @param material_reader Reader that provides access to the material data of the input file 
    * @param unit_handler Instance to provide (non-)dimensionalization of values 
    * @return Full initialized material pairing
    * 
    * @note during the initialization checks are done to read only variables that are required with the given 
    *       compile time settings 
    */
   MaterialPairing InitializeMaterialPairing( std::vector<unsigned int> const& material_indices, 
                                              MaterialReader const& material_reader,
                                              UnitHandler const& unit_handler ) {

      //**************************************************************************************************************************************
      // 1. Read all data from input file 
      // Declare all variables that can be filled (use default values)
      double surface_tension_coefficient_fixed_value = -1.0;

      // Surface tension coefficient only if needed
      if constexpr( CC::CapillaryForcesActive() ) {
         // read model data if required  
         if constexpr( CC::SurfaceTensionCoefficientModelActive() ) {
            throw std::runtime_error( "Currently no surface tension coefficient model is implemented!" );
         } else {// otherwise read the fixed value 
            surface_tension_coefficient_fixed_value = material_reader.ReadFixedValue( material_indices, MaterialProperty::SurfaceTensionCoefficient );
         }
      } 

      //**************************************************************************************************************************************
      // 2. Create final data (eos + modles) and log information
      // logger
      LogWriter & logger = LogWriter::Instance();  
      logger.LogMessage( " " );
      logger.LogMessage( "Material pairing " + std::to_string( material_indices[0] ) + " <-> " + std::to_string( material_indices[1] ) + ":" );

      // Material pairing properties 
      // Log properties only if positive values are specified (they exist)
      std::string tmp_string;
      // Bulk viscosity
      tmp_string = StringOperations::Indent( 2 ) + "Surface tension coefficient: ";
      if( surface_tension_coefficient_fixed_value < 0.0 ) {
         logger.LogMessage( tmp_string + "Not Used" );
      } else {
         logger.LogMessage( tmp_string + StringOperations::ToScientificNotationString( surface_tension_coefficient_fixed_value, 9 ) );
      }
      
      //**************************************************************************************************************************************
      // 3. return the fully initialied material pairing
      return MaterialPairing( surface_tension_coefficient_fixed_value, unit_handler );
   }

} // namespace initialization