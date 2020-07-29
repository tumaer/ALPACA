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
#include "initialization/initialization_modular_algorithm_assembler.h"

#include "user_specifications/compile_time_constants.h"
#include "utilities/string_operations.h"

namespace Initialization {

   /**
    * @brief Computes the gravity source term.
    * @param source_term_reader The reader providing access to all source term information of the input file.
    * param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return gravity in x-, y- and z-direction.
    */
   std::array<double, 3> GetGravity( SourceTermReader const& source_term_reader, UnitHandler const& unit_handler ) {
      // Initialize the gravity with zero
      std::array<double, 3> gravity = { 0.0, 0.0, 0.0 };

      // Only proceed if gravity is activated
      if constexpr( CC::GravityIsActive() ) {
         gravity[0] = unit_handler.NonDimensionalizeValue( source_term_reader.ReadGravity( Direction::X ), { UnitType::Length }, { UnitType::Time, UnitType::Time } );
         // For two and three dimensions
         if constexpr( CC::DIM() != Dimension::One ) {
            gravity[1] = unit_handler.NonDimensionalizeValue( source_term_reader.ReadGravity( Direction::Y ), { UnitType::Length }, { UnitType::Time, UnitType::Time } );
         }
         // Only for three dimensions
         if constexpr( CC::DIM() == Dimension::Three ) {
            gravity[2] = unit_handler.NonDimensionalizeValue( source_term_reader.ReadGravity( Direction::Z ), { UnitType::Length }, { UnitType::Time, UnitType::Time } );
         }
      }

      return gravity;
   }

   /**
    * @brief Get a vector containing all levels of the current simulation.
    * @param maximum_level Maximum level used for the simulation.
    * @return vector with all levels.
    */
   std::vector<unsigned int> GetAllLevels( unsigned int const maximum_level ) {
      std::vector<unsigned int> all_levels( maximum_level + 1 );//Level zero need to be counted as well
      std::iota( all_levels.begin(), all_levels.end(), 0 );
      return all_levels;
   }

   /**
    * @brief Initializes the complete modular algorithm assembler class with the given input classes.
    * @param input_reader Reader that provides access to the full data of the input file.
    * @param topology_manager Class providing global (on all ranks) node information.
    * @param tree Tree class providing local (on current rank) node information.
    * @param communication_manager Calls providing communication handling between different ranks.
    * @param multiresolution Instance to provide mutliresolution computations for remeshing and so fourth.
    * @param material_manager Instance providing initialized material data.
    * @param input_output_manager Instance to provide restart handling and output writing.
    * @param initial_condition Instance for handling initial conditions.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return The fully initialized Modular algorithm assembler class.
    */
   ModularAlgorithmAssembler InitializeModularAlgorithmAssembler( InputReader const& input_reader,
                                                                  TopologyManager& topology_manager,
                                                                  Tree& tree,
                                                                  CommunicationManager& communication_manager,
                                                                  HaloManager& halo_manager,
                                                                  Multiresolution const& multiresolution,
                                                                  MaterialManager const& material_manager,
                                                                  InputOutputManager& input_output_manager,
                                                                  InitialCondition const& initial_condition,
                                                                  UnitHandler const& unit_handler ) {

      // Get data that is logged
      double const start_time = unit_handler.NonDimensionalizeValue( input_reader.GetTimeControlReader().ReadStartTime(), UnitType::Time );
      double const end_time   = unit_handler.NonDimensionalizeValue( input_reader.GetTimeControlReader().ReadEndTime(), UnitType::Time );
      double const cfl_number = input_reader.GetTimeControlReader().ReadCFLNumber();

      // Log data
      LogWriter& logger = LogWriter::Instance();
      logger.LogMessage( " " );
      logger.LogMessage( "Time control data: " );
      logger.LogMessage( StringOperations::Indent( 2 ) + "Start time: " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( start_time, UnitType::Time ), 9 ) );
      logger.LogMessage( StringOperations::Indent( 2 ) + "End time  : " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( end_time, UnitType::Time ), 9 ) );
      logger.LogMessage( StringOperations::Indent( 2 ) + "CFL number: " + StringOperations::ToScientificNotationString( cfl_number, 9 ) );
      logger.LogMessage( " " );
      // Compute the cell size on maximum level
      unsigned int const maximum_level = topology_manager.GetMaximumLevel();
      // The MAA required the dimensionless cell size -> No dimensionalization required
      double const cell_size_on_maximum_level = tree.GetNodeSizeOnLevelZero() / double( CC::ICX() ) / double( 1 << maximum_level );

      // initialize the algorithm assembler
      return ModularAlgorithmAssembler( start_time,
                                        end_time,
                                        cfl_number,
                                        GetGravity( input_reader.GetSourceTermReader(), unit_handler ),
                                        GetAllLevels( maximum_level ),
                                        cell_size_on_maximum_level,
                                        unit_handler,
                                        initial_condition,
                                        tree,
                                        topology_manager,
                                        halo_manager,
                                        communication_manager,
                                        multiresolution,
                                        material_manager,
                                        input_output_manager );
   }
}// namespace Initialization