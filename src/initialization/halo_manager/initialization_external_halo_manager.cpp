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
#include "initialization/halo_manager/initialization_external_halo_manager.h"

#include "prime_states/prime_state_handler_setup.h"

namespace Initialization {

   /**
    * @brief Converts the user given (reduced) set of primestates into a full set of conservatives for each material.
    * @param fixed_input_prime_states Input (reduced) set of prime states.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @param material_manager material_manager Instance providing initialized material data.
    * @return Vector with the fixed conservatives for all materials.
    */
   std::vector<std::array<double, MF::ANOE()>> ConvertInputToConservatives( std::array<double, MF::ANOP()> const& fixed_input_prime_states,
                                                                            UnitHandler const& unit_handler,
                                                                            MaterialManager const& material_manager ) {

      // Declare the prime state handler that is used in the simulation and initialize it
      using PrimeStateHandlerConcretization = PrimeStateHandlerSetup::Concretize<prime_state_handler>::type;
      PrimeStateHandlerConcretization const prime_state_handler( material_manager );
      // Get all materials contained in the simulation
      std::vector<MaterialName> const material_names( material_manager.GetMaterialNames() );

      // non-dimensionalize prime states with temprorary object
      std::array<double, MF::ANOP()> fixed_prime_states;
      for( PrimeState const p : MF::ASOP() ) {
         fixed_prime_states[PTI( p )] = unit_handler.NonDimensionalizeValue( fixed_input_prime_states[PTI( p )], MF::FieldUnit( p ) );
      }

      // Vector that is returned holding all conservative fixed values for each material
      std::vector<std::array<double, MF::ANOE()>> fixed_conservatives;
      fixed_conservatives.resize( material_names.size() );
      //initialize fixed conservatives at boundary with 0
      for( unsigned int mat_index = 0; mat_index < material_names.size(); ++mat_index ) {
         for( unsigned int conservative_index = 0; conservative_index < MF::ANOE(); ++conservative_index ) {
            fixed_conservatives[mat_index][conservative_index] = 0.0;
         }
      }
      // Compute conservatives from inputfile prime states
      for( unsigned int mat_index = 0; mat_index < material_names.size(); ++mat_index ) {
         prime_state_handler.ConvertPrimeStatesToConservatives( material_names[mat_index], fixed_prime_states, fixed_conservatives[mat_index] );
      }

      return fixed_conservatives;
   }

   /**
    * @brief Converts the full set of fixed conservatives computed from the ( reduced ) prime states into the full set of prime states.
    * @return Vector with full set of fixed prime states for all materials.
    */
   std::vector<std::array<double, MF::ANOP()>> ConvertConservativesToPrimeStates( std::vector<std::array<double, MF::ANOE()>> const& fixed_conservatives,
                                                                                  MaterialManager const& material_manager ) {

      // Declare the prime state handler that is used in the simulation and initialize it
      using PrimeStateHandlerConcretization = PrimeStateHandlerSetup::Concretize<prime_state_handler>::type;
      PrimeStateHandlerConcretization const prime_state_handler( material_manager );

      // Get all materials contained in the simulation
      std::vector<MaterialName> const material_names( material_manager.GetMaterialNames() );
      // Vector that is returned (reserved enough memory)
      std::vector<std::array<double, MF::ANOP()>> fixed_prime_states;
      fixed_prime_states.resize( material_names.size() );
      //initialize fixed prime_states at boundary with 0
      for( unsigned int material_index = 0; material_index < material_names.size(); ++material_index ) {
         for( unsigned int primestate_index = 0; primestate_index < MF::ANOP(); ++primestate_index ) {
            fixed_prime_states[material_index][primestate_index] = 0.0;
         }
      }

      // compute active set of prime states from inputfile prime states
      for( unsigned int mat_index = 0; mat_index < material_names.size(); ++mat_index ) {
         prime_state_handler.ConvertConservativesToPrimeStates( material_names[mat_index], fixed_conservatives[mat_index], fixed_prime_states[mat_index] );
      }

      return fixed_prime_states;
   }

   /**
   * @brief Initializes the material boundary conditions with the given input data. Periodic Boundary conditions are not created.
   * @param bc_reader Reader to provide access to the boundary conditions in the input data.
   * @param unit_handler Instance to provide (non-)dimensionalization of values.
   * @param material_manager material_manager Instance providing initialized material data.
   * @return Array with pointers to the base class for all material boundary conditions.
   */
   std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6> InitializeMaterialBoundaryConditions( BoundaryConditionReader const& bc_reader,
                                                                                                         UnitHandler const& unit_handler,
                                                                                                         MaterialManager const& material_manager ) {
      // Declare vector that is returned
      std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6> material_boundary_conditions;
      // Initialize with nullptr
      for( unsigned int i = 0; i < 6; i++ ) {
         material_boundary_conditions[i] = nullptr;
      }

      // Create logger for input logging
      LogWriter & logger = LogWriter::Instance();

      logger.LogMessage( " " );
      logger.LogMessage( StringOperations::Indent( 2 ) + "Materials:" );

      // x- direction
      material_boundary_conditions[LTI( BoundaryLocation::East )] = CreateMaterialBoundary<BoundaryLocation::East>( bc_reader, unit_handler, material_manager );
      material_boundary_conditions[LTI( BoundaryLocation::West )] = CreateMaterialBoundary<BoundaryLocation::West>( bc_reader, unit_handler, material_manager );
      // y-direction
      if constexpr( CC::DIM() != Dimension::One ) {
         material_boundary_conditions[LTI( BoundaryLocation::South )] = CreateMaterialBoundary<BoundaryLocation::South>( bc_reader, unit_handler, material_manager );
         material_boundary_conditions[LTI( BoundaryLocation::North )] = CreateMaterialBoundary<BoundaryLocation::North>( bc_reader, unit_handler, material_manager );
      }

      // z-direction
      if constexpr( CC::DIM() == Dimension::Three ) {
         material_boundary_conditions[LTI( BoundaryLocation::Top )]    = CreateMaterialBoundary<BoundaryLocation::Top>( bc_reader, unit_handler, material_manager );
         material_boundary_conditions[LTI( BoundaryLocation::Bottom )] = CreateMaterialBoundary<BoundaryLocation::Bottom>( bc_reader, unit_handler, material_manager );
      }

      return material_boundary_conditions;
   }

   /**
   * @brief Initializes the levelset boundary conditions with the given input data. Periodic Boundary conditions are not created.
   * @param bc_reader Reader to provide access to the boundary conditions in the input data.
   * @return Array with pointers to the base class for all levelset boundary conditions.
   *
   * @note The levelset boundary conditions are always required, even if single fluid simulations are done. This is due to the Interface tag updates
   *       that are done for single and multi-material simulations.
   */
   std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> InitializeLevelsetBoundaryConditions( BoundaryConditionReader const& bc_reader ) {
      // Declare vector that is returned
      std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> levelset_boundary_conditions;
      // Initialize with nullptr
      for( unsigned int i = 0; i < 6; i++ ) {
         levelset_boundary_conditions[i] = nullptr;
      }

      // Create logger for input logging
      LogWriter & logger = LogWriter::Instance();

      logger.LogMessage( " " );
      logger.LogMessage( StringOperations::Indent( 2 ) +  "Levelset:" );

      // x- direction
      levelset_boundary_conditions[LTI( BoundaryLocation::East )] = CreateLevelsetBoundary<BoundaryLocation::East>( bc_reader );
      levelset_boundary_conditions[LTI( BoundaryLocation::West )] = CreateLevelsetBoundary<BoundaryLocation::West>( bc_reader );

      // y-direction
      if constexpr( CC::DIM() != Dimension::One ) {
         levelset_boundary_conditions[LTI( BoundaryLocation::South )] = CreateLevelsetBoundary<BoundaryLocation::South>( bc_reader );
         levelset_boundary_conditions[LTI( BoundaryLocation::North )] = CreateLevelsetBoundary<BoundaryLocation::North>( bc_reader );
      }

      // z-direction
      if constexpr( CC::DIM() == Dimension::Three ) {
         levelset_boundary_conditions[LTI( BoundaryLocation::Top )]    = CreateLevelsetBoundary<BoundaryLocation::Top>( bc_reader );
         levelset_boundary_conditions[LTI( BoundaryLocation::Bottom )] = CreateLevelsetBoundary<BoundaryLocation::Bottom>( bc_reader );
      }

      return levelset_boundary_conditions;
   }

   /**
    * @brief Initializes the complete external halo manager class with the given input reader and other classes.
    * @param input_reader Reader that provides access to the full data of the input file.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @param material_manager material_manager Instance providing initialized material data.
    * @return The fully initialized ExternalHaloManager class.
    */
   ExternalHaloManager InitializeExternalHaloManager( InputReader const& input_reader,
                                                      UnitHandler const& unit_handler,
                                                      MaterialManager const& material_manager ) {

      // Create logger for input logging
      LogWriter & logger = LogWriter::Instance();

      logger.LogMessage( " " );
      logger.LogMessage( "External boundary conditions:" );

      // First create and then move to provide proper logging of data (correct order)
      std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6> material_boundaries(
            InitializeMaterialBoundaryConditions( input_reader.GetBoundaryConditionReader(), unit_handler, material_manager ) );

      std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> levelset_boundaries(
            InitializeLevelsetBoundaryConditions( input_reader.GetBoundaryConditionReader() ) );

      logger.LogMessage( " " );

      // Note: Here, ownership transfer takes place. First initialized pointers are now nullptrs
      return ExternalHaloManager( std::move( material_boundaries ), std::move( levelset_boundaries ) );
   }


   /**
    * @brief Provides a proper string for the data of the fixed value boundary condition.
    * @param indent Indention width used for the logging string.
    * @param fixed_prime_states Values of the fixed prime states for all materials.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return The full log string.
    */
   std::string LogFixedValueData( unsigned int const indent,
                                  std::vector<std::array<double, MF::ANOP()>> const fixed_prime_states,
                                  UnitHandler const& unit_handler ) {

      std::string tmp_string;

      // Loop through all active prime states to find the variable with longest name
      std::size_t maximum_size = 0;
      std::vector<std::string> prime_names( MF::ASOP().size() );
      for( PrimeState const& prime : MF::ASOP() ) {
         // Get the defined input name
         std::string prime_string = std::string( MF::InputName( prime ) );
         // If name does not exist use specific prefix with the index number
         if( prime_string.empty() ) prime_string = "Primestate_" + std::to_string( PTI( prime ) + 1 );
         // add the name to the vector
         prime_names[PTI( prime )] = prime_string;
         // Determine maximum size
         maximum_size = prime_string.size() > maximum_size ? prime_string.size() : maximum_size;
      }

      // Loop through all materials and prime states to print the data
      for( size_t mat_index = 0; mat_index < fixed_prime_states.size(); mat_index++ ) {
         // Write the string of the material (+1 is used since the materials start with 1 in input file but with zero in the vector)
         tmp_string += StringOperations::Indent( indent + 2 ) + "Material " + std::to_string( mat_index + 1 ) + "\n";
         // Add all dimensionalized values of the primestates
         for( PrimeState const& prime : MF::ASOP() ) {
            tmp_string +=   StringOperations::Indent( indent + 4 ) + prime_names[PTI( prime )] + std::string( maximum_size - prime_names[PTI( prime )].size(), ' ' ) + ": "
                          + StringOperations::ToScientificNotationString(
                              unit_handler.DimensionalizeValue( fixed_prime_states[mat_index][PTI( prime )], MF::FieldUnit( prime ) ), 9, true ) + "\n";
         }
         tmp_string += " \n";
      }
      return tmp_string;
   }
}