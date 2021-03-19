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
#include "instantiation/instantiation_initial_condition.h"

#include <stdexcept>
#include "utilities/string_operations.h"

#include "initial_condition/levelset_initializer.h"
#include "initial_condition/functional_levelset_initializer.h"
#include "initial_condition/stl_levelset_initializer.h"
#include "initial_condition/parametric_levelset_initializer.h"

namespace Instantiation {

   /**
    * @brief Reads the initial material condition expression string for all materials.
    * @param initial_condition_reader Instance that provides access to the initial condition data in the input file.
    * @param number_of_materials Number of materials present in the current simulation.
    * @param prime_state_variable_names All names of the prime states that are read (contains empty strings for prime states that should not be read).
    * @return Initial condition expression strings for all materials.
    */
   std::vector<std::string> GetMaterialInitialConditions( InitialConditionReader const& initial_condition_reader,
                                                          unsigned int const number_of_materials,
                                                          std::vector<std::string> const& prime_state_variable_names ) {
      // logger for warning logging
      LogWriter& logger = LogWriter::Instance();

      // Vector that is returned
      std::vector<std::string> material_initial_conditions( number_of_materials );
      // Loop through all materials (start at one )
      for( unsigned int material_index = 0; material_index < number_of_materials; material_index++ ) {
         // material index +1 is used since in input file the indices start at 1 and internally at 0
         material_initial_conditions[material_index] = initial_condition_reader.ReadMaterialInitialConditions( material_index + 1 );
         // Check if all non-empty prime state variables are given. Otherwise throw a warning
         for( auto const& variable_name : prime_state_variable_names ) {
            if( !variable_name.empty() ) {
               if( material_initial_conditions[material_index].find( variable_name ) == std::string::npos ) {
                  logger.LogMessage( " " );
                  logger.LogMessage( "Warning!! The variable name '" + variable_name + "' is not contained in the initial condition string for material " + std::to_string( material_index + 1 ) + "!" );
                  logger.LogMessage( " " );
               }
            }
         }
      }

      return material_initial_conditions;
   }

   /**
    * @brief Gives the variables for the parametric initializer.
    * @param initial_condition_reader Instance to read the initial condition data.
    * @param levelset_index The levelset index for which the variables should be read.
    * @return The array of parametric variables.
    */
   std::array<ParametricVariable, 2> CreateParametricVariables( InitialConditionReader const& initial_condition_reader, unsigned int const levelset_index ) {
      // Read all variables
      std::vector<ParametricVariable> const parametric_variables( initial_condition_reader.ReadParametricLevelsetInitializerVariables( levelset_index ) );
      // Assign the correct number of variables depending on the given input an used dimensions
      if( parametric_variables.size() >= 2 && CC::DIM() == Dimension::Three ) {
         return { parametric_variables[0], parametric_variables[1] };
      } else if( parametric_variables.size() >= 1 ) {
         return { parametric_variables[0], ParametricVariable() };
      } else {
         return { ParametricVariable(), ParametricVariable() };
      }
   }

   /**
    * @brief Returns a pointer to the levelset initializer obtained from the given input data.
    * @param initial_condition_reader Instance to read the initial condition data.
    * @param levelset_index The levelset index for which the variables should be read.
    * @param material_names Names of the materials in the simulation.
    * @param node_size_on_level_zero_ Dimensional size of each node on level zero.
    * @param maximum_level Maximum level of the simulation.
    * @return Unique pointer to the const base class of the levelset initializer.
    */
   std::unique_ptr<LevelsetInitializer> InstantiateLevelsetInitializer( InitialConditionReader const& initial_condition_reader,
                                                                        unsigned int const levelset_index,
                                                                        std::vector<MaterialName> const& material_names,
                                                                        double const node_size_on_level_zero_,
                                                                        unsigned int const maximum_level ) {

      // If the number of materials is 1 (single material) do not require reading and assign a constant expression (cannot be empty. UserExpression throws error otherwise.)
      if( material_names.size() == 1 ) {
         return std::make_unique<FunctionalLevelsetInitializer>( "phi := 1.0;", std::vector<std::array<double, 6>>(), material_names, node_size_on_level_zero_, maximum_level );
      }
      // read the type and input of the levelset initializer
      LevelsetInitializerType const initializer_type{ initial_condition_reader.ReadLevelsetInitializerType( levelset_index, LevelsetInitializerType::Functional ) };
      std::string const initializer_input{ initial_condition_reader.ReadLevelsetInitializerInput( levelset_index ) };
      // read the bounding boxex
      std::vector<std::array<double, 6>> bounding_boxes{ initial_condition_reader.ReadLevelsetInitializerBoundingBoxes( levelset_index ) };

      // logger
      LogWriter& logger = LogWriter::Instance();
      logger.LogMessage( " " );
      // input data
      logger.LogMessage( "Levelset initializer: " );

      // switch between different levelset initializer types to call specific constructor
      switch( initializer_type ) {
         case LevelsetInitializerType::Functional: {
            // 1. Create, 2. Log, 3. Return levelset_initializer
            std::unique_ptr<LevelsetInitializer> levelset_initializer( std::make_unique<FunctionalLevelsetInitializer>( initializer_input, bounding_boxes, material_names, node_size_on_level_zero_, maximum_level ) );
            logger.LogMessage( levelset_initializer->GetLogData( 2 ) );
            logger.LogMessage( " " );
            return levelset_initializer;
         }
         case LevelsetInitializerType::STL: {
            // 1. Create, 2. Log, 3. Return levelset_initializer
            std::unique_ptr<LevelsetInitializer> levelset_initializer( std::make_unique<StlLevelsetInitializer>( StringOperations::RemoveSpaces( initializer_input ), bounding_boxes, material_names, node_size_on_level_zero_, maximum_level ) );
            logger.LogMessage( levelset_initializer->GetLogData( 2 ) );
            logger.LogMessage( " " );
            return levelset_initializer;
         }
         case LevelsetInitializerType::Parametric: {
            // For this type of initializer we need to read additional data
            std::array<ParametricVariable, 2> const parametric_variables( CreateParametricVariables( initial_condition_reader, levelset_index ) );
            std::array<double, 3> const reference_point( initial_condition_reader.ReadParametricLevelsetInitializerReferencePoint( levelset_index ) );
            // 1. Create, 2. Log, 3. Return levelset_initializer
            std::unique_ptr<LevelsetInitializer> levelset_initializer( std::make_unique<ParametricLevelsetInitializer>( initializer_input, parametric_variables, reference_point,
                                                                                                                        bounding_boxes, material_names, node_size_on_level_zero_, maximum_level ) );
            logger.LogMessage( levelset_initializer->GetLogData( 2 ) );
            logger.LogMessage( " " );
            return levelset_initializer;
         }
         default: {
            throw std::logic_error( "The chosen levelset initialization method has not been implemented!" );
         }
      }
   }

   /**
    * @brief Instantiates the complete initial condition class with the given input classes.
    * @param input_reader Reader that provides access to the full data of the input file.
    * @param topology_manager Class providing global (on all ranks) node information.
    * @param tree Tree class providing local (on current rank) node information.
    * @param material_manager Instance providing instantiated material data.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return The fully instantiated InititalCondition class.
    */
   std::unique_ptr<InitialCondition> InstantiateInitialCondition( InputReader const& input_reader,
                                                                  TopologyManager const& topology_manager,
                                                                  Tree const& tree,
                                                                  MaterialManager const& material_manager,
                                                                  UnitHandler const& unit_handler ) {

      // Create the prime states variable names
      // Get all variables that used in the initial condition string (defined in field details)
      std::vector<std::string> prime_state_variable_names( MF::ANOP() );
      for( PrimeState const p : MF::ASOP() ) {
         // No check is done here if variable should be used for input. Later 0 is assigned to those
         prime_state_variable_names[PTI( p )] = StringOperations::RemoveSpaces( std::string( MF::InputName( p ) ) );
      }

      //For more than one initializer this must be a loop and the resulting initializer a vector
      unsigned int const levelset_index = 0;

      std::unique_ptr<LevelsetInitializer> levelset_initializer( InstantiateLevelsetInitializer( input_reader.GetInitialConditionReader(),
                                                                                                 levelset_index + 1,
                                                                                                 material_manager.GetMaterialNames(),
                                                                                                 unit_handler.DimensionalizeValue( tree.GetNodeSizeOnLevelZero(), UnitType::Length ),
                                                                                                 topology_manager.GetMaximumLevel() ) );

      // Initialize the prime state handler
      std::vector<std::string> const material_initial_conditions{ GetMaterialInitialConditions( input_reader.GetInitialConditionReader(), material_manager.GetNumberOfMaterials(), prime_state_variable_names ) };
      std::unique_ptr<PrimeStateInitializer const> prime_state_initializer{ std::make_unique<PrimeStateInitializer const>( material_initial_conditions,
                                                                                                                           prime_state_variable_names,
                                                                                                                           unit_handler.DimensionalizeValue( tree.GetNodeSizeOnLevelZero(), UnitType::Length ),
                                                                                                                           unit_handler ) };
      return std::make_unique<InitialCondition>( std::move( prime_state_initializer ),
                                                 std::move( levelset_initializer ) );
   }

}// namespace Instantiation
