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
#include "instantiation/materials/instantiation_material_pairing.h"

#include "user_specifications/compile_time_constants.h"
#include "materials/material_property_definitions.h"
#include "materials/material_pairing_property_models/surface_tension_coefficient_models/constant_surface_tension_coefficient_model.h"

namespace Instantiation {

   /**
    * @brief Gives the model for the material pairing property surface tension coefficient for the given input data.
    * @param model_name Name of the surface tension coefficient model.
    * @param model_data Data of the surface tension coefficient model.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return pointer to the const base class of all material parameter models.
    */
   std::unique_ptr<InterfaceParameterModel const> InstantiateSurfaceTensionCoefficientModel( MaterialPropertyModelName const model_name,
                                                                                             std::unordered_map<std::string, double> const& model_data,
                                                                                             UnitHandler const& unit_handler ) {
      // logger
      LogWriter& logger = LogWriter::Instance();

      // switch between different surface tension coefficient models
      switch( model_name ) {
         case MaterialPropertyModelName::NotUsed: {
            return nullptr;
         }
         case MaterialPropertyModelName::SurfaceTensionCoefficientConstant: {
            // 1. Create, 2. Log, 3. Return model
            std::unique_ptr<ConstantSurfaceTensionCoefficientModel const> model( std::make_unique<ConstantSurfaceTensionCoefficientModel const>( model_data, unit_handler ) );
            logger.LogLinebreakMessage( model->GetLogData( 6, unit_handler ) );
            return model;
         }
         default: {
            throw std::logic_error( "This surface tension coefficient model is not implemented yet!" );
         }
      }
   }

   /**
    * @brief Instantiates a complete material pairing with the given input reader.
    * @param material_indices Indices of the materials that form the pairing, which should be read from the input file.
    * @param material_reader Reader that provides access to the material data of the input file.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return Full initialized material pairing.
    *
    * @note during the initialization checks are done to read only variables that are required with the given
    *       compile time settings.
    */
   MaterialPairing InstantiateMaterialPairing( std::vector<unsigned int> const& material_indices,
                                               MaterialReader const& material_reader,
                                               UnitHandler const& unit_handler ) {

      // Read all data from input file
      // Declare all variables that can be filled (use default values)
      double surface_tension_coefficient_fixed_value                                   = -1.0;
      MaterialPropertyModelName surface_tension_coefficient_model_name                 = MaterialPropertyModelName::NotUsed;
      std::unordered_map<std::string, double> surface_tension_coefficient_model_data   = {};
      std::unique_ptr<InterfaceParameterModel const> surface_tension_coefficient_model = nullptr;

      // Surface tension coefficient only if needed
      if constexpr( CC::CapillaryForcesActive() ) {
         // read model data if required
         if constexpr( CC::SurfaceTensionCoefficientModelActive() ) {
            // read the model data
            surface_tension_coefficient_model_name = material_reader.ReadModelName( material_indices, MaterialProperty::SurfaceTensionCoefficient );
            surface_tension_coefficient_model_data = material_reader.ReadModelData( material_indices, MaterialProperty::SurfaceTensionCoefficient );
         } else {// otherwise read the fixed value
            surface_tension_coefficient_fixed_value = material_reader.ReadFixedValue( material_indices, MaterialProperty::SurfaceTensionCoefficient );
         }
      }

      // Create final data (eos + modles) and log information
      // logger
      LogWriter& logger = LogWriter::Instance();
      logger.LogMessage( " " );
      logger.LogMessage( "Material pairing " + std::to_string( material_indices[0] ) + " <-> " + std::to_string( material_indices[1] ) + ":" );

      // Material pairing properties
      // Log properties only if positive values are specified (they exist)
      std::string tmp_string;
      // Bulk viscosity
      tmp_string = StringOperations::Indent( 2 ) + "Surface tension coefficient: ";
      if constexpr( CC::CapillaryForcesActive() ) {
         if( surface_tension_coefficient_model_name != MaterialPropertyModelName::NotUsed ) {
            // create model (logs data itself)
            logger.LogMessage( tmp_string );
            surface_tension_coefficient_model = InstantiateSurfaceTensionCoefficientModel( surface_tension_coefficient_model_name, surface_tension_coefficient_model_data, unit_handler );
         } else {
            logger.LogMessage( tmp_string + StringOperations::ToScientificNotationString( surface_tension_coefficient_fixed_value, 9 ) );
         }
      } else {
         logger.LogMessage( tmp_string + "Not Used" );
      }

      // return the fully initialied material pairing
      return MaterialPairing( surface_tension_coefficient_fixed_value, std::move( surface_tension_coefficient_model ), unit_handler );
   }

}// namespace Instantiation
