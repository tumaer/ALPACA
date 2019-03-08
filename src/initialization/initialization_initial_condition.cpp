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
#include "initialization/initialization_initial_condition.h"

#include <stdexcept>
#include "utilities/string_operations.h"
namespace Initialization {

   /**
    * @brief Reads the initial material condition expression string for all materials
    * @param initial_condition_reader Instance that provides access to the initial condition data in the input file
    * @param number_of_materials Number of materials present in the current simulation
    * @param variable_names_prime_states All names of the prime states that are read (contains empty strings for prime states that should not be read)
    * @return Initial condition expression strings for all materials
    */
   std::vector<std::string> GetMaterialInitialConditions( InitialConditionReader const& initial_condition_reader,
                                                          unsigned int const number_of_materials,
                                                          std::vector<std::string> const& variable_names_prime_states ) {
      // logger for warning logging
      LogWriter & logger = LogWriter::Instance();

      // Vector that is returned
      std::vector<std::string> material_initial_conditions( number_of_materials );
      // Loop through all materials (start at one )
      for( unsigned int material_index = 0; material_index < number_of_materials; material_index++ ) {
         // material index +1 is used since in input file the indices start at 1 and internally at 0
         material_initial_conditions[material_index] = initial_condition_reader.ReadMaterialInitialConditions( material_index + 1 );
         // Check if all non-empty prime state variables are given. Otherwise throw a warning
         for( auto const& variable_name : variable_names_prime_states ) {
            if( !variable_name.empty() ) {
               if( material_initial_conditions[material_index].find( variable_name ) == std::string::npos ) {
                  logger.LogMessage( " " );
                  logger.LogMessage(   "Warning!! The variable name '" + variable_name + "' is not contained in the initial condition string for material "
                                     + std::to_string( material_index + 1 ) + "!" );
                  logger.LogMessage( " " );
               }
            }
         }
      }

      return material_initial_conditions;
   }

   /**
    * @brief Reads the initial levelset condition expressions of all interfaces
    * @param initial_condition_reader Instance that provides access to the initial condition data in the input file
    * @param number_of_materials Number of materials present in the current simulation
    * @param variable_name_levelset Name of the levelset variable that is used
    * @return Initial condition expression strings for all levelsets
    */
   std::vector<std::string> GetLevelsetInitialConditions( InitialConditionReader const& initial_condition_reader,
                                                          unsigned int number_of_materials,
                                                          std::string const& variable_name_levelset ) {
      // If the number of materials is 1 (single material) do not require reading
      if( number_of_materials == 1 ) {
         return {};
      }

      // logger for warning logging
      LogWriter & logger = LogWriter::Instance();

      // Vector that is returned (one element less than number of materials since all levelset initial conditions are referred to material 1)
      std::vector<std::string> levelset_initial_conditions( number_of_materials - 1 );
      // Loop through all levelsets (start at one since input file starts at 1). Only N-1 interfaces are present 1 <-> 2, 1 <-> 3, 1 <-> 4, ...
      for( unsigned int levelset_index = 0; levelset_index < number_of_materials - 1; levelset_index++ ) {
         levelset_initial_conditions[levelset_index] = initial_condition_reader.ReadLevelsetInitialConditions( levelset_index + 1 );
         // Check if levelset variables is given in the string. Otherwise throw a warning
         if( levelset_initial_conditions[levelset_index].find( variable_name_levelset ) == std::string::npos ) {
            logger.LogMessage( " " );
            logger.LogMessage(   "Warning!! The levelset variable name '" + variable_name_levelset + "' is not contained in the initial condition string for levelset "
                               + std::to_string( levelset_index + 1 ) + "!" );
            logger.LogMessage( " " );
         }
      }

      return levelset_initial_conditions;
   }

   /**
    * @brief Initializes the complete initial condition class with the given input classes
    * @param input_reader Reader that provides access to the full data of the input file
    * @param topology_manager Class providing global (on all ranks) node information
    * @param tree Tree class providing local (on current rank) node information
    * @param material_manager Instance providing initialized material data
    * @param unit_handler Instance to provide (non-)dimensionalization of values
    * @return The fully initialized InititalCondition class
    */
   InitialCondition InitializeInitialCondition( InputReader const& input_reader,
                                                TopologyManager const& topology_manager,
                                                Tree const& tree,
                                                MaterialManager const& material_manager,
                                                UnitHandler const& unit_handler ) {

      // Create the prime states variable names
      // Get all variables that used in the initial condition string (defined in field details)
      std::vector<std::string> variable_names_prime_states( MF::ANOP() );
      for( PrimeState const p : MF::ASOP() ) {
         // No check is done here if variable should be used for input. Later 0 is assigned to those
         variable_names_prime_states[PTI( p )] = StringOperations::RemoveSpaces( std::string( MF::InputName( p ) ) );
      }

      // Create the levelset variable in case multi fluid simulations are used
      std::string variable_name_levelset = "";
      if( material_manager.GetNumberOfMaterials() > 1 ) {
         variable_name_levelset = StringOperations::RemoveSpaces( std::string( IF::InputName( InterfaceDescription::Levelset ) ) );
         // Check if empty and throw error in that case
         if( variable_name_levelset.empty() ) {
            throw std::logic_error( "The input name of the levelset field must be given for multi material simulations!" );
         }
      }

      return InitialCondition( GetMaterialInitialConditions( input_reader.GetInitialConditionReader(), material_manager.GetNumberOfMaterials(), variable_names_prime_states ),
                               GetLevelsetInitialConditions( input_reader.GetInitialConditionReader(), material_manager.GetNumberOfMaterials(), variable_name_levelset ),
                               material_manager.GetMaterialNames(),
                               variable_names_prime_states,
                               variable_name_levelset,
                               unit_handler.DimensionalizeValue( tree.GetNodeSizeOnLevelZero(), UnitType::Length ),
                               topology_manager.GetMaximumLevel(),
                               unit_handler );
   }

} // namespace Initialization