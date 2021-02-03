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
#ifndef INITIALIZATION_EXTERNAL_HALO_MANAGER_H
#define INITIALIZATION_EXTERNAL_HALO_MANAGER_H

#include <vector>
#include <array>

#include "input_output/input_reader.h"
#include "unit_handler.h"
#include "materials/material_manager.h"

#include "boundary_condition/external_halo_manager.h"
#include "boundary_condition/boundary_specifications.h"
#include "boundary_condition/material_boundary_condition.h"
#include "boundary_condition/levelset_boundary_condition.h"
#include "boundary_condition/symmetry_boundary_condition.h"
#include "boundary_condition/zero_gradient_boundary_condition.h"
#include "boundary_condition/wall_boundary_condition.h"
#include "boundary_condition/fixed_value_boundary_condition.h"

/**
 * @brief Defines all instantiation functions required for the external halo manager.
 */
namespace Instantiation {

   // Initialization of the full external halo manager
   ExternalHaloManager InstantiateExternalHaloManager( InputReader const& input_reader,
                                                       UnitHandler const& unit_handler,
                                                       MaterialManager const& material_manager );

   // Initialization of the full set of material boundary conditions
   std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6> InstantiateMaterialBoundaryConditions( BoundaryConditionReader const& bc_reader,
                                                                                                          UnitHandler const& unit_handler,
                                                                                                          MaterialManager const& material_manager );
   // Initialization of the full set of levelset boundary conditions
   std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> InstantiateLevelsetBoundaryConditions( BoundaryConditionReader const& bc_reader );

   // Converts the (reduced) set of prime states into a full set of conservatives
   std::vector<std::array<double, MF::ANOE()>> ConvertInputToConservatives( std::array<double, MF::ANOP()> const& input_fixed_prime_states,
                                                                            UnitHandler const& unit_handler,
                                                                            MaterialManager const& material_manager );

   // Converts the full set of conservatives into a full set of prime states
   std::vector<std::array<double, MF::ANOP()>> ConvertConservativesToPrimeStates( std::vector<std::array<double, MF::ANOE()>> const& fixed_conservatives,
                                                                                  MaterialManager const& material_manager );

   // Returns a string for logging fixed value data
   std::string LogFixedValueData( unsigned int const indent,
                                  std::vector<std::array<double, MF::ANOP()>> const fixed_prime_states,
                                  UnitHandler const& unit_handler );

   /**
    * @brief Builds a material boundary condition of appropiate type.
    * @param bc_reader The boundary condition reader providing acces to boundary information in the input file.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @param material_manager Instance providing initialized material data.
    *
    * @tparam LOC Template parameter for the boundary location.
    */
   template<BoundaryLocation LOC>
   inline std::unique_ptr<MaterialBoundaryCondition const> CreateMaterialBoundary( BoundaryConditionReader const& bc_reader,
                                                                                   UnitHandler const& unit_handler,
                                                                                   MaterialManager const& material_manager ) {
      // Create logger for input logging
      LogWriter& logger = LogWriter::Instance();

      // Obtain the correct material boundary type from the reader
      MaterialBoundaryType const material_boundary_type( bc_reader.ReadMaterialBoundaryType( LOC ) );

      //Switch with return do not need break.
      switch( material_boundary_type ) {
         case MaterialBoundaryType::ZeroGradient: {
            // 1. Create the boundary condition, 2. Log its information, 3. return it
            std::unique_ptr<ZeroGradientBoundaryCondition<LOC> const> boundary_condition( std::make_unique<ZeroGradientBoundaryCondition<LOC> const>() );

            logger.LogMessage( StringOperations::Indent( 4 ) + BoundaryLocationToString( LOC, true ) + std::string( BoundaryLocationToString( BoundaryLocation::Bottom, true ).size() - BoundaryLocationToString( LOC, true ).size(), ' ' ) + ": Zero-gradient" );

            return boundary_condition;
         }
         case MaterialBoundaryType::Symmetry: {
            // 1. Create the boundary condition, 2. Log its information, 3. return it
            std::unique_ptr<SymmetryBoundaryCondition<LOC> const> boundary_condition( std::make_unique<SymmetryBoundaryCondition<LOC> const>() );

            logger.LogMessage( StringOperations::Indent( 4 ) + BoundaryLocationToString( LOC, true ) + std::string( BoundaryLocationToString( BoundaryLocation::Bottom, true ).size() - BoundaryLocationToString( LOC, true ).size(), ' ' ) + ": Symmetry" );

            return boundary_condition;
         }
         case MaterialBoundaryType::FixedValue: {
            // NOTE: Approach:
            //        1. take the (potentially reduced) set of prime states defined by the user in the input file
            //        2. Convert the reduced set into a full set of conservatives
            //        3. Revert the full set of conservatives into a full set of prime states

            // Obtain the (reduced) set of fixed prime states from the input file
            std::array<double, MF::ANOP()> const input_fixed_prime_states( bc_reader.ReadMaterialFixedValueBoundaryConditions( LOC ) );
            // Compute the fixed value conservatives from the reduced set
            std::vector<std::array<double, MF::ANOE()>> const fixed_value_conservatives(
                  ConvertInputToConservatives( input_fixed_prime_states, unit_handler, material_manager ) );
            // Compute the full set of fixed value prime states
            std::vector<std::array<double, MF::ANOP()>> const fixed_value_prime_states(
                  ConvertConservativesToPrimeStates( fixed_value_conservatives, material_manager ) );

            // 1. Create the fixed boundary condition, 2. Log its information, 3. return it
            std::unique_ptr<FixedValueBoundaryCondition<LOC> const> boundary_condition(
                  std::make_unique<FixedValueBoundaryCondition<LOC> const>( fixed_value_conservatives, fixed_value_prime_states ) );

            logger.LogMessage( StringOperations::Indent( 4 ) + BoundaryLocationToString( LOC, true ) + std::string( BoundaryLocationToString( BoundaryLocation::Bottom, true ).size() - BoundaryLocationToString( LOC, true ).size(), ' ' ) + ": Fixed value" );
            logger.LogMessage( " " );
            logger.LogLinebreakMessage( LogFixedValueData( 4, fixed_value_prime_states, unit_handler ) );

            return boundary_condition;
         }
         case MaterialBoundaryType::Wall: {
            // 1. Create the boundary condition, 2. Log its information, 3. return it
            std::unique_ptr<WallBoundaryCondition<LOC> const> boundary_condition( std::make_unique<WallBoundaryCondition<LOC> const>() );

            logger.LogMessage( StringOperations::Indent( 4 ) + BoundaryLocationToString( LOC, true ) + std::string( BoundaryLocationToString( BoundaryLocation::Bottom, true ).size() - BoundaryLocationToString( LOC, true ).size(), ' ' ) + ": Wall" );

            return boundary_condition;
         }
         case MaterialBoundaryType::Periodic: {
            // Not handled here. Simply log it
            logger.LogMessage( StringOperations::Indent( 4 ) + BoundaryLocationToString( LOC, true ) + std::string( BoundaryLocationToString( BoundaryLocation::Bottom, true ).size() - BoundaryLocationToString( LOC, true ).size(), ' ' ) + ": Periodic" );

            return nullptr;
         }
         case MaterialBoundaryType::Internal: {
            throw std::logic_error( "Internals Halos are not managed by the ExternalHaloManager!" );
         }
         default:
            throw std::invalid_argument( "No such Material Boundary Condition exists" );
      }
   }

   /**
    * @brief Builds a levelset boundary condition of appropiate type.
    * @param bc_reader The boundary condition reader providing acces to boundary information in the input file.
    *
    * @tparam LOC Template parameter for the boundary location.
    */
   template<BoundaryLocation LOC>
   inline std::unique_ptr<LevelsetBoundaryCondition const> CreateLevelsetBoundary( BoundaryConditionReader const& bc_reader ) {

      // Create logger for input logging
      LogWriter& logger = LogWriter::Instance();

      // Obtain the correct levelset boundary type from the reader
      LevelSetBoundaryType const levelset_boundary_type( bc_reader.ReadLevelsetBoundaryType( LOC ) );

      //Switch with return do not need break.
      switch( levelset_boundary_type ) {
         case LevelSetBoundaryType::ZeroGradient: {
            // 1. Create the boundary condition, 2. Log its information, 3. return it
            std::unique_ptr<ZeroGradientBoundaryCondition<LOC> const> boundary_condition( std::make_unique<ZeroGradientBoundaryCondition<LOC> const>() );

            logger.LogMessage( StringOperations::Indent( 4 ) + BoundaryLocationToString( LOC, true ) + std::string( BoundaryLocationToString( BoundaryLocation::Bottom, true ).size() - BoundaryLocationToString( LOC, true ).size(), ' ' ) + ": Zero-gradient" );

            return boundary_condition;
         }
         case LevelSetBoundaryType::Symmetry: {
            // 1. Create the boundary condition, 2. Log its information, 3. return it
            std::unique_ptr<SymmetryBoundaryCondition<LOC> const> boundary_condition( std::make_unique<SymmetryBoundaryCondition<LOC> const>() );

            logger.LogMessage( StringOperations::Indent( 4 ) + BoundaryLocationToString( LOC, true ) + std::string( BoundaryLocationToString( BoundaryLocation::Bottom, true ).size() - BoundaryLocationToString( LOC, true ).size(), ' ' ) + ": Symmetry" );

            return boundary_condition;
         }
         case LevelSetBoundaryType::Periodic: {
            // Not handled here. Simply log it
            logger.LogMessage( StringOperations::Indent( 4 ) + BoundaryLocationToString( LOC, true ) + std::string( BoundaryLocationToString( BoundaryLocation::Bottom, true ).size() - BoundaryLocationToString( LOC, true ).size(), ' ' ) + ": Periodic" );

            return nullptr;
         }
         case LevelSetBoundaryType::Internal: {
            throw std::logic_error( "Internals Halos are not managed by the ExternalHaloManager!" );
         }
         default: {
            throw std::invalid_argument( "No such Levelset Boundary Condition exists" );
         }
      }
   }
}// namespace Instantiation

#endif// INITIALIZATION_EXTERNAL_HALO_MANAGER_H
