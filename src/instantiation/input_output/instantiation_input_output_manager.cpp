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
#include "instantiation/input_output/instantiation_input_output_manager.h"

#include "input_output/utilities/file_utilities.h"
#include "utilities/string_operations.h"

namespace Instantiation {

   /**
    * @brief Computes the output time stamps to write a certain output.
    * @param output_reader Reader that provide access to the output specific input data.
    * @param time_control_reader Reader that provide access to the time control specific input data.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @param output_type Type of the output that should be used (standard, debug, interface).
    * @return Vector with all output time stamps for the given output type.
    */
   std::vector<double> ComputeOutputTimes( OutputReader const& output_reader,
                                           TimeControlReader const& time_control_reader,
                                           UnitHandler const& unit_handler,
                                           OutputType const output_type ) {

      // Declare the vector to be returned
      std::vector<double> output_times;
      // Read the time type (Off, Interval or Stamps)
      OutputTimesType const times_type( output_reader.ReadOutputTimesType( output_type ) );
      // Read the start and end time and dimensionalize them (checks on consistency are done in reader)
      double const start_time( unit_handler.NonDimensionalizeValue( time_control_reader.ReadStartTime(), UnitType::Time ) );
      double const end_time( unit_handler.NonDimensionalizeValue( time_control_reader.ReadEndTime(), UnitType::Time ) );

      // Off
      if( times_type == OutputTimesType::Off ) {
         // Return empty vector immediately to pass the checks done afterwards
         return output_times;
      }

      // Interval times type
      if( times_type == OutputTimesType::Interval || times_type == OutputTimesType::IntervalStamps ) {
         // Read the desired interval of times (check on consistency done in reader)
         double const interval( unit_handler.NonDimensionalizeValue( output_reader.ReadOutputInterval( output_type ), UnitType::Time ) );
         // fill output_times with timestamps generated from given interval (ceil ensures that at least the number of time stamps is 1)
         size_t number_of_timesteps = size_t( std::ceil( std::nextafter( ( end_time - start_time ) / interval, 0 ) ) );

         // Remove the virtually generated step for intervals larger than start_time - end_time to exclude end_time (will be added in the end exactly)
         if( number_of_timesteps != 0 ) {
            number_of_timesteps -= 1;
         }

         // Resize vector to its actual number of steps
         output_times.resize( number_of_timesteps );
         // Fill the vector with all values from start_time in the given interval
         double time = start_time;
         std::generate( output_times.begin(), output_times.end(), [&time, &interval]() { return ( time += interval ); } );
      }

      // Timestamps times type
      if( times_type == OutputTimesType::Stamps || times_type == OutputTimesType::IntervalStamps ) {
         // Possible that both options are used. Therefore, first create local and then add them to the final timestamps vector
         std::vector<double> output_times_tmp;

         // Read the time stamps from the input file
         output_times_tmp = output_reader.ReadOutputTimeStamps( output_type );
         // Dimensionalize each time stamp (must be done before comparison to start and end time)
         std::for_each( output_times_tmp.begin(), output_times_tmp.end(),
                        [&unit_handler]( double& timestamp ) { timestamp = unit_handler.NonDimensionalizeValue( timestamp, UnitType::Time ); } );
         // Remove values smaller equal than start time
         output_times_tmp.erase( std::remove_if( output_times_tmp.begin(), output_times_tmp.end(),
                                                 [&start_time]( double const& timestamp ) { return timestamp <= start_time; } ),
                                 output_times_tmp.end() );
         // Remove values larger equal than end time
         output_times_tmp.erase( std::remove_if( output_times_tmp.begin(), output_times_tmp.end(),
                                                 [&end_time]( double const& timestamp ) { return timestamp >= end_time; } ),
                                 output_times_tmp.end() );

         // Append all time stamps to the global vector and sort them
         output_times.reserve( output_times.size() + output_times_tmp.size() );
         output_times.insert( output_times.end(), output_times_tmp.begin(), output_times_tmp.end() );
         std::sort( output_times.begin(), output_times.end() );
      }

      // now add the exact end_time as final timestamp
      output_times.push_back( end_time );

      // Final checks (could be removed since all requirements are already considered)
      if( std::any_of( output_times.begin(), output_times.end(), [&start_time]( double const& timestamp ) { return timestamp < start_time; } ) ) {
         throw std::invalid_argument( "Output timestamps for " + OutputTypeToString( output_type ) + " output must be larger than start_time" );
      }
      if( std::any_of( output_times.begin(), output_times.end(), [&end_time]( double const& timestamp ) { return timestamp > end_time; } ) ) {
         throw std::invalid_argument( "Output timestamps for " + OutputTypeToString( output_type ) + " output must be smaller than end_time" );
      }
      if( !std::is_sorted( output_times.begin(), output_times.end() ) ) {
         throw std::invalid_argument( "Output timestamps for " + OutputTypeToString( output_type ) + " output must be ascending" );
      }

      return output_times;
   }

   /**
    * @brief Computes the restart time stamps to write a certain restart snapshot file.
    * @param restart_reader Reader that provide access to the restart specific input data.
    * @param time_control_reader Reader that provide access to the time control specific input data.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return Vector with all snaposhot time stamps.
    */
   std::vector<double> ComputeSnapshotTimes( RestartReader const& restart_reader, TimeControlReader const& time_control_reader, UnitHandler const& unit_handler ) {
      // Declare vector to be returned
      std::vector<double> snapshot_times;
      // read the start and end time
      double const start_time = unit_handler.NonDimensionalizeValue( time_control_reader.ReadStartTime(), UnitType::Time );
      double const end_time   = unit_handler.NonDimensionalizeValue( time_control_reader.ReadEndTime(), UnitType::Time );

      // Read the time type
      SnapshotTimesType const times_type( restart_reader.ReadSnapshotTimesType() );

      // Compute times stamps only if flag is set
      if( times_type == SnapshotTimesType::Stamps || times_type == SnapshotTimesType::IntervalStamps ) {
         snapshot_times = restart_reader.ReadSnapshotTimeStamps();
         std::for_each( snapshot_times.begin(), snapshot_times.end(), [&unit_handler]( double& timestamp ) { timestamp = unit_handler.NonDimensionalizeValue( timestamp, UnitType::Time ); } );
         // Remove all time stamps that are smaller than the start time
         snapshot_times.erase( std::remove_if( snapshot_times.begin(), snapshot_times.end(),
                                               [&start_time]( double const& timestamp ) { return timestamp < start_time; } ),
                               snapshot_times.end() );
         // Remove all time stamps that are larger than the end time
         snapshot_times.erase( std::remove_if( snapshot_times.begin(), snapshot_times.end(),
                                               [&end_time]( double const& timestamp ) { return timestamp > end_time; } ),
                               snapshot_times.end() );
         // Add the exact end time explicitly
         snapshot_times.push_back( end_time );

         // final checks
         if( std::any_of( snapshot_times.begin(), snapshot_times.end(), [&start_time]( double const& timestamp ) { return timestamp < start_time; } ) ) {
            throw std::invalid_argument( "Snapshot timestamps must be larger than start_time" );
         }
         if( std::any_of( snapshot_times.begin(), snapshot_times.end(), [&end_time]( double const& timestamp ) { return timestamp > end_time; } ) ) {
            throw std::invalid_argument( "Snapshot timestamps must be smaller than end_time" );
         }
         if( !std::is_sorted( snapshot_times.begin(), snapshot_times.end() ) ) {
            throw std::invalid_argument( "Snapshot timestamps must be ascending" );
         }
      }

      return snapshot_times;
   }

   /**
    * @brief Instantiates the complete input output manager class with the given input reader and other input classes.
    * @param input_reader Reader that provides access to the full data of the input file.
    * @param topology_manager Class providing global (on all ranks) node information.
    * @param tree Tree class providing local (on current rank) node information.
    * @param material_manager Instance providing initialized material data.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return The fully instantiated MaterialManager class.
    */
   InputOutputManager InstantiateInputOutputManager( InputReader const& input_reader,
                                                     OutputWriter const& output_writer,
                                                     RestartManager const& restart_manager,
                                                     UnitHandler const& unit_handler ) {

      // Get the required readers
      TimeControlReader const& time_control_reader( input_reader.GetTimeControlReader() );
      RestartReader const& restart_reader( input_reader.GetRestartReader() );
      OutputReader const& output_reader( input_reader.GetOutputReader() );

      // Read all data that is required
      // Input
      InputType const input_type   = input_reader.GetInputType();
      std::string const input_file = input_reader.GetInputFile();
      // Output
      std::string output_folder_name                        = FileUtilities::AddUnusedNumberToPath( FileUtilities::RemoveFilePath( FileUtilities::RemoveFileExtension( input_file ) ) );
      std::vector<double> const standard_output_timestamps  = ComputeOutputTimes( output_reader, time_control_reader, unit_handler, OutputType::Standard );
      std::vector<double> const interface_output_timestamps = ComputeOutputTimes( output_reader, time_control_reader, unit_handler, OutputType::Interface );
      double const time_naming_factor                       = output_reader.ReadTimeNamingFactor();
      // Restart
      RestoreMode const restore_mode              = restart_reader.ReadRestoreMode();
      std::string const restart_file              = restore_mode != RestoreMode::Off ? restart_reader.ReadRestoreFilename() : "";
      SnapshotTimesType const snapshot_times_type = restart_reader.ReadSnapshotTimesType();
      std::vector<double> snapshot_timestamps     = ComputeSnapshotTimes( restart_reader, time_control_reader, unit_handler );
      unsigned int const snapshot_interval        = ( snapshot_times_type == SnapshotTimesType::Interval || snapshot_times_type == SnapshotTimesType::IntervalStamps ) ?
                                                   restart_reader.ReadSnapshotInterval() :
                                                   0;
      unsigned int const snapshots_to_keep = ( snapshot_times_type == SnapshotTimesType::Interval || snapshot_times_type == SnapshotTimesType::IntervalStamps ) ?
                                                   restart_reader.ReadSnapshotIntervalsToKeep() :
                                                   0;

      // logging
      LogWriter& logger = LogWriter::Instance();
      logger.LogMessage( " " );
      // input data
      logger.LogMessage( "File/Folder information: " );
      logger.LogMessage( StringOperations::Indent( 2 ) + "Input type        : " + InputTypeToString( input_type ) );
      logger.LogMessage( StringOperations::Indent( 2 ) + "Simulation Name   : " + FileUtilities::RemoveFilePath( FileUtilities::RemoveFileExtension( input_file ) ) );
      logger.LogMessage( StringOperations::Indent( 2 ) + "Output Folder     : " + output_folder_name );
      logger.LogMessage( StringOperations::Indent( 2 ) + "Time naming factor: " + StringOperations::ToScientificNotationString( time_naming_factor, 9 ) );
      logger.LogMessage( " " );
      // Output and restart
      // Header for time stamps
      if( !standard_output_timestamps.empty() || !interface_output_timestamps.empty() || !snapshot_timestamps.empty() ) {
         logger.LogMessage( StringOperations::Indent( 37 ) + "First" + StringOperations::Indent( 10 ) + "Last" + StringOperations::Indent( 6 ) + "Total" );
      }
      // output time stamps
      if( !standard_output_timestamps.empty() || !interface_output_timestamps.empty() ) {
         if( !standard_output_timestamps.empty() ) {
            logger.LogMessage( "Standard output time stamps :    " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( standard_output_timestamps.front(), UnitType::Time ), 6 ) + "  " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( standard_output_timestamps.back(), UnitType::Time ), 6 ) + "    " + std::to_string( standard_output_timestamps.size() ) );
         }
         if( !interface_output_timestamps.empty() ) {
            logger.LogMessage( "Interface output time stamps:    " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( interface_output_timestamps.front(), UnitType::Time ), 6 ) + "  " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( interface_output_timestamps.back(), UnitType::Time ), 6 ) + "    " + std::to_string( interface_output_timestamps.size() ) );
         }
      }
      // restart time stamps
      logger.LogMessage( " " );
      if( !snapshot_timestamps.empty() ) {
         logger.LogMessage( "Snapshot time stamps        :    " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( snapshot_timestamps.front(), UnitType::Time ), 6 ) + "  " + StringOperations::ToScientificNotationString( unit_handler.DimensionalizeValue( snapshot_timestamps.back(), UnitType::Time ), 6 ) + "    " + std::to_string( snapshot_timestamps.size() ) );
      }
      if( snapshot_interval > 0 ) {
         logger.LogMessage( "Restart snapshot interval   : " + std::to_string( snapshot_interval ) );
         logger.LogMessage( "Kept interval snapshots     : " + std::to_string( snapshots_to_keep ) );
      }

      if( snapshot_interval == 0 && snapshot_timestamps.empty() ) {
         logger.LogMessage( "Restart snapshots           : Disabled" );
      }

      if( standard_output_timestamps.empty() ) {
         logger.LogMessage( "Standard output files       : Disabled" );
      } else if( interface_output_timestamps.empty() ) {
         logger.LogMessage( "Interface output files      : Disabled" );
      }

      logger.LogMessage( " " );

      // Instantiate the input output manager
      return InputOutputManager( input_file,
                                 output_folder_name,
                                 unit_handler,
                                 output_writer,
                                 restart_manager,
                                 time_naming_factor,
                                 standard_output_timestamps,
                                 interface_output_timestamps,
                                 restore_mode,
                                 restart_file,
                                 snapshot_timestamps,
                                 snapshot_interval,
                                 snapshots_to_keep );
   }

}// namespace Instantiation
