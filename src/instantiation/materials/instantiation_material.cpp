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
#include "instantiation/materials/instantiation_material.h"

#include "materials/equation_of_state.h"
#include "user_specifications/compile_time_constants.h"

#include "materials/equations_of_state/stiffened_gas.h"
#include "materials/equations_of_state/stiffened_gas_safe.h"
#include "materials/equations_of_state/stiffened_gas_complete_safe.h"
#include "materials/equations_of_state/waterlike_fluid.h"
#include "materials/equations_of_state/noble_abel_stiffened_gas.h"
#include "materials/equations_of_state/isentropic.h"

#include "materials/material_property_models/shear_viscosity_models/constant_shear_viscosity_model.h"
#include "materials/material_property_models/shear_viscosity_models/shear_rate_models/power_law_shear_viscosity_model.h"
#include "materials/material_property_models/shear_viscosity_models/shear_rate_models/cross_shear_viscosity_model.h"
#include "materials/material_property_models/shear_viscosity_models/shear_rate_models/carreau_yasuda_shear_viscosity_model.h"
#include "materials/material_property_models/shear_viscosity_models/temperature_models/sutherland_shear_viscosity_model.h"
#include "materials/material_property_models/thermal_conductivity_models/constant_thermal_conductivity_model.h"

namespace Instantiation {

   /**
    * @brief Returns a pointer to the equation of state obtained from the given input data.
    * @param eos_name Name of the equation of state.
    * @param eos_data All parameters specified for the equation of state from the input file.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return pointer to the const base class of all equations of state.
    */
   std::unique_ptr<EquationOfState const> InstantiateEquationOfState( EquationOfStateName const eos_name,
                                                                      std::unordered_map<std::string, double> const& eos_data,
                                                                      UnitHandler const& unit_handler ) {
      // logger
      LogWriter& logger = LogWriter::Instance();

      // switch between different eos to call specific constructors
      switch( eos_name ) {
         case EquationOfStateName::StiffenedGas: {
            // 1. Create, 2. Log, 3. Return eos
            std::unique_ptr<StiffenedGas const> eos( std::make_unique<StiffenedGas const>( eos_data, unit_handler ) );
            logger.LogLinebreakMessage( eos->GetLogData( 4, unit_handler ) );
            return eos;
         }
         case EquationOfStateName::StiffenedGasSafe: {
            // 1. Create, 2. Log, 3. Return eos
            std::unique_ptr<StiffenedGasSafe const> eos( std::make_unique<StiffenedGasSafe const>( eos_data, unit_handler ) );
            logger.LogLinebreakMessage( eos->GetLogData( 4, unit_handler ) );
            return eos;
         }
         case EquationOfStateName::StiffenedGasCompleteSafe: {
            // 1. Create, 2. Log, 3. Return eos
            std::unique_ptr<StiffenedGasCompleteSafe const> eos( std::make_unique<StiffenedGasCompleteSafe const>( eos_data, unit_handler ) );
            logger.LogLinebreakMessage( eos->GetLogData( 4, unit_handler ) );
            return eos;
         }
         case EquationOfStateName::WaterlikeFluid: {
            // 1. Create, 2. Log, 3. Return eos
            std::unique_ptr<WaterlikeFluid const> eos( std::make_unique<WaterlikeFluid const>( eos_data, unit_handler ) );
            logger.LogLinebreakMessage( eos->GetLogData( 4, unit_handler ) );
            return eos;
         }
         case EquationOfStateName::NobleAbelStiffenedGas: {
            // sanity check if Gruneisen flag ist
            if constexpr( !CC::GruneisenDensityDependent() ) {
               throw std::runtime_error( "To use NobleAbelStiffenedGas you need to activate CC::GruneisenDensityDependent" );
            }
            // 1. Create, 2. Log, 3. Return eos
            std::unique_ptr<NobleAbelStiffenedGas const> eos( std::make_unique<NobleAbelStiffenedGas const>( eos_data, unit_handler ) );
            logger.LogLinebreakMessage( eos->GetLogData( 4, unit_handler ) );
            return eos;
         }
         case EquationOfStateName::Isentropic: {
            // 1. Create, 2. Log, 3. Return eos
            std::unique_ptr<Isentropic const> eos( std::make_unique<Isentropic const>( eos_data, unit_handler ) );
            logger.LogLinebreakMessage( eos->GetLogData( 4, unit_handler ) );
            return eos;
         }
         default: {
            throw std::logic_error( "The equation of state has not been implemented!" );
         }
      }
   }

   /**
    * @brief Gives the model for the material property shear viscosity for the given input data.
    * @param model_name Name of the shear viscosity model.
    * @param model_data Data of the shear viscosity model.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return pointer to the const base class of all material parameter models.
    */
   std::unique_ptr<MaterialParameterModel const> InstantiateShearViscosityModel( MaterialPropertyModelName const model_name,
                                                                                 std::unordered_map<std::string, double> const& model_data,
                                                                                 UnitHandler const& unit_handler ) {
      // logger
      LogWriter& logger = LogWriter::Instance();

      // switch between different shear viscosity models
      switch( model_name ) {
         case MaterialPropertyModelName::NotUsed: {
            return nullptr;
         }
         case MaterialPropertyModelName::ShearViscosityConstant: {
            // 1. Create, 2. Log, 3. Return model
            std::unique_ptr<ConstantShearViscosityModel const> model( std::make_unique<ConstantShearViscosityModel const>( model_data, unit_handler ) );
            logger.LogLinebreakMessage( model->GetLogData( 6, unit_handler ) );
            return model;
         }
         case MaterialPropertyModelName::ShearViscosityPowerLaw: {
            // 1. Create, 2. Log, 3. Return model
            std::unique_ptr<PowerLawShearViscosityModel const> model( std::make_unique<PowerLawShearViscosityModel const>( model_data, unit_handler ) );
            logger.LogLinebreakMessage( model->GetLogData( 6, unit_handler ) );
            return model;
         }
         case MaterialPropertyModelName::ShearViscosityCarreauYasuda: {
            // 1. Create, 2. Log, 3. Return model
            std::unique_ptr<CarreauYasudaShearViscosityModel const> model( std::make_unique<CarreauYasudaShearViscosityModel const>( model_data, unit_handler ) );
            logger.LogLinebreakMessage( model->GetLogData( 6, unit_handler ) );
            return model;
         }
         case MaterialPropertyModelName::ShearViscosityCross: {
            // 1. Create, 2. Log, 3. Return model
            std::unique_ptr<CrossShearViscosityModel const> model( std::make_unique<CrossShearViscosityModel const>( model_data, unit_handler ) );
            logger.LogLinebreakMessage( model->GetLogData( 6, unit_handler ) );
            return model;
         }
         case MaterialPropertyModelName::ShearViscositySutherland: {
            // 1. Create, 2. Log, 3. Return model
            std::unique_ptr<SutherlandShearViscosityModel const> model( std::make_unique<SutherlandShearViscosityModel const>( model_data, unit_handler ) );
            logger.LogLinebreakMessage( model->GetLogData( 6, unit_handler ) );
            return model;
         }
         default: {
            throw std::logic_error( "This shear viscosity model is not implemented yet!" );
         }
      }
   }

   /**
    * @brief Gives the model for the material property thermal conductivity for the given input data.
    * @param model_name Name of the thermal conductivity model.
    * @param model_data Data of the thermal conductivity model.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return pointer to the const base class of all material parameter models.
    */
   std::unique_ptr<MaterialParameterModel const> InstantiateThermalConductivityModel( MaterialPropertyModelName const model_name,
                                                                                      std::unordered_map<std::string, double> const& model_data,
                                                                                      UnitHandler const& unit_handler ) {
      // logger
      LogWriter& logger = LogWriter::Instance();

      // switch between different conductivity models
      switch( model_name ) {
         case MaterialPropertyModelName::NotUsed: {
            return nullptr;
         }
         case MaterialPropertyModelName::ThermalConductivityConstant: {
            // 1. Create, 2. Log, 3. Return model
            std::unique_ptr<ConstantThermalConductivityModel const> model( std::make_unique<ConstantThermalConductivityModel const>( model_data, unit_handler ) );
            logger.LogLinebreakMessage( model->GetLogData( 6, unit_handler ) );
            return model;
         }
         default: {
            throw std::logic_error( "This thermal conductivity model is not implemented yet!" );
         }
      }
   }

   /**
    * @brief Instantiates the material type and makes consistency checks.
    * @param material_index Index of the material to be read from the input file.
    * @param material_reader Reader that provides access to the material data of the input file.
    * @return The checked material type.
    */
   MaterialType InstantiateMaterialType( unsigned int const material_index, MaterialReader const& material_reader ) {
      // Read the material type (default to Fluid if not given)
      MaterialType const material_type( material_reader.ReadMaterialType( material_index, MaterialType::Fluid ) );

      // Make consistency checks, that the type corresponds to a type that is allowed with current user specifications
      if constexpr( !CC::SolidBoundaryActive() ) {
         if( material_type == MaterialType::SolidBoundary ) {
            throw std::runtime_error( "To use solid boundary materials active the compile time constants." );
         }
      }

      return material_type;
   }
   /* @brief Instantiates a complete material with the given input reader.* @param material_index Index of the material to be read from the input file.* @param material_reader Reader that provides access to the material data of the input file.* @param unit_handler Instance to provide( non - ) dimensionalization of values.* @ return Full initialized material.*
         * @note during the initialization checks are done to read only variables that are required with the given
         * compile time settings.
         */
   std::tuple<MaterialType, Material> InstantiateMaterial( unsigned int const material_index,
                                                           MaterialReader const& material_reader,
                                                           UnitHandler const& unit_handler ) {
      // Read all data from input file
      // Read the material type (default to Fluid if not given)
      MaterialType const material_type( InstantiateMaterialType( material_index + 1, material_reader ) );

      // Always read the data for the equation of state
      EquationOfStateName const eos_name( material_reader.ReadEquationOfStateName( material_index + 1 ) );
      // data map cannot be const since specific heat is added below in case for NobelAbelStiffened Gas
      std::unordered_map<std::string, double> eos_data( material_reader.ReadEquationOfStateData( material_index + 1 ) );

      // Declare all variables that can be filled (use default values)
      double specific_heat_capacity_fixed_value                                = -1.0;
      double thermal_conductivity_fixed_value                                  = -1.0;
      double bulk_viscosity_fixed_value                                        = -1.0;
      double shear_viscosity_fixed_value                                       = -1.0;
      MaterialPropertyModelName shear_viscosity_model_name                     = MaterialPropertyModelName::NotUsed;
      std::unordered_map<std::string, double> shear_viscosity_model_data       = {};
      std::unique_ptr<MaterialParameterModel const> shear_viscosity_model      = nullptr;
      MaterialPropertyModelName thermal_conductivity_model_name                = MaterialPropertyModelName::NotUsed;
      std::unordered_map<std::string, double> thermal_conductivity_model_data  = {};
      std::unique_ptr<MaterialParameterModel const> thermal_conductivity_model = nullptr;

      // Declare the vector that is used for reading the correct values
      // material index +1 is used since in input file the indices start at 1 and internally at 0
      std::vector<unsigned int> const material_indices = { material_index + 1 };

      // Read viscosities only if required
      if constexpr( CC::ViscosityIsActive() ) {
         // always read the fixed value bulk viscosity
         bulk_viscosity_fixed_value = material_reader.ReadFixedValue( material_indices, MaterialProperty::BulkViscosity );
         // shear viscosity
         // read model data if required
         if constexpr( CC::ShearViscosityModelActive() ) {
            // read the model data
            shear_viscosity_model_name = material_reader.ReadModelName( material_indices, MaterialProperty::ShearViscosity );
            shear_viscosity_model_data = material_reader.ReadModelData( material_indices, MaterialProperty::ShearViscosity );
         } else {// otherwise read the fixed value
            shear_viscosity_fixed_value = material_reader.ReadFixedValue( material_indices, MaterialProperty::ShearViscosity );
         }
      }

      // thermal conductivity only if needed
      if constexpr( CC::HeatConductionActive() ) {
         // thermal conductivity model data
         if constexpr( CC::ThermalConductivityModelActive() ) {
            // read the model data
            thermal_conductivity_model_name = material_reader.ReadModelName( material_indices, MaterialProperty::ThermalConductivity );
            thermal_conductivity_model_data = material_reader.ReadModelData( material_indices, MaterialProperty::ThermalConductivity );
         } else {
            // thermal conductivity fixed value
            thermal_conductivity_fixed_value = material_reader.ReadFixedValue( material_indices, MaterialProperty::ThermalConductivity );
         }
      }

      // specific heat (property for heat conduction)
      if constexpr( CC::HeatConductionActive() ) {
         specific_heat_capacity_fixed_value = material_reader.ReadFixedValue( material_indices, MaterialProperty::SpecificHeatCapacity );
      }

      // specific heat (property for equation of state). Second read cycle to not fill default specific heat fixed value.
      if( eos_name == EquationOfStateName::NobleAbelStiffenedGas ) {
         std::string const specific_heat_name( MaterialPropertyToString( MaterialProperty::SpecificHeatCapacity, false ) );
         eos_data[specific_heat_name] = material_reader.ReadFixedValue( material_indices, MaterialProperty::SpecificHeatCapacity );
      }

      // Create final data (type + eos + models) and log information
      // logger
      LogWriter& logger = LogWriter::Instance();
      logger.LogMessage( " " );
      logger.LogMessage( "Material " + std::to_string( material_index + 1 ) + ":" );

      // Material type
      logger.LogMessage( std::string( 60, '-' ) );
      logger.LogMessage( StringOperations::Indent( 2 ) + "Material type: " + MaterialTypeToString( material_type ) );

      // Equation of state
      logger.LogMessage( std::string( 60, '-' ) );
      logger.LogMessage( StringOperations::Indent( 2 ) + "Equation of state:" );

      // create eos
      std::unique_ptr<EquationOfState const> eos( InstantiateEquationOfState( eos_name, eos_data, unit_handler ) );

      // Material properties
      logger.LogMessage( std::string( 60, '-' ) );
      logger.LogMessage( StringOperations::Indent( 2 ) + "Properties:" );
      // Log properties only if positive values are specified (they exist)
      std::string tmp_string;
      // Bulk viscosity
      tmp_string = StringOperations::Indent( 4 ) + "Bulk viscosity        : ";
      if constexpr( CC::ViscosityIsActive() ) {
         logger.LogMessage( tmp_string + StringOperations::ToScientificNotationString( bulk_viscosity_fixed_value, 9 ) );
      } else {
         logger.LogMessage( tmp_string + "Not Used" );
      }
      // Shear viscosity
      tmp_string = StringOperations::Indent( 4 ) + "Shear viscosity       : ";
      if constexpr( CC::ViscosityIsActive() ) {
         if( shear_viscosity_model_name != MaterialPropertyModelName::NotUsed ) {
            // create model (logs data itself)
            logger.LogMessage( tmp_string );
            shear_viscosity_model = InstantiateShearViscosityModel( shear_viscosity_model_name, shear_viscosity_model_data, unit_handler );
         } else {
            logger.LogMessage( tmp_string + StringOperations::ToScientificNotationString( shear_viscosity_fixed_value, 9 ) );
         }
      } else {
         logger.LogMessage( tmp_string + "Not Used" );
      }
      // Thermal conductivity
      tmp_string = StringOperations::Indent( 4 ) + "Thermal conductivity  : ";
      if constexpr( CC::HeatConductionActive() ) {
         if( thermal_conductivity_model_name != MaterialPropertyModelName::NotUsed ) {
            // create model (logs data itself)
            logger.LogMessage( tmp_string );
            thermal_conductivity_model = InstantiateThermalConductivityModel( thermal_conductivity_model_name, thermal_conductivity_model_data, unit_handler );
         } else {
            logger.LogMessage( tmp_string + StringOperations::ToScientificNotationString( thermal_conductivity_fixed_value, 9 ) );
         }
      } else {
         logger.LogMessage( tmp_string + "Not Used" );
      }
      // Specific heat capacity
      tmp_string = StringOperations::Indent( 4 ) + "Specific heat capacity: ";
      if constexpr( CC::HeatConductionActive() ) {
         logger.LogMessage( tmp_string + StringOperations::ToScientificNotationString( specific_heat_capacity_fixed_value, 9 ) );
      } else {
         logger.LogMessage( tmp_string + "Not Used" );
      }

      logger.LogMessage( std::string( 60, '-' ) );

      // Return the final fully instantiated material ( move operations to transfer pointer ownership (deleted move constructor))
      return std::make_tuple( material_type, Material( std::move( eos ),
                                                       bulk_viscosity_fixed_value,
                                                       shear_viscosity_fixed_value,
                                                       thermal_conductivity_fixed_value,
                                                       specific_heat_capacity_fixed_value,
                                                       std::move( shear_viscosity_model ),
                                                       std::move( thermal_conductivity_model ),
                                                       unit_handler ) );
   }
}// namespace Instantiation
