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
#include "input_output_manager.h"

#include <cstdio> // needed for file deletion
#include <fstream>
#include <unistd.h>
#include <algorithm>

#include "utilities/file_operations.h"

#include "user_specifications/debug_and_profile_setup.h"
#include "input_output/output_writer/output_definitions.h"
#include "input_output/output_writer.h"
#include "user_specifications/compile_time_constants.h"
#include "communication/mpi_utilities.h"

// Helper function for better readability in the code below
namespace {
   inline void RemoveTimeStamps( std::vector<double> & time_stamps, double const maxmimum_value ) {
      time_stamps.erase( time_stamps.begin(), std::find_if_not(time_stamps.begin(), time_stamps.end(), [&maxmimum_value]( double const timestamp )
                                                               { return timestamp <= maxmimum_value; } ) );
   }
}

/**
 * @brief Default constructor using the inputs from the operating system.
 * @param input_file Input file name with the full path
 * @param output_folder Name of the output folder where the data is written into
 * @param unit_handler Instance to provide (non-)dimensionalization of values
 * @param output_writer Full initialized output writer class
 * @param restart_manager Full initialized restar_manager class
 * @param time_naming_factor Factor that is used for naming the ouput files
 * @param standard_output_timestamps Timestamps when output is desired for the standard or debug output
 * @param interface_output_timestamps Timestamps when output is desired for the interface output
 * @param restore_mode Restore mode (Off, soft, forced)
 * @param restore_filename Filename that is used for restoring the simulation if mode is enabled
 * @param restart_snapshot_timestampsTimestamps when restart files are desired
 * @param restart_snapshot_interval Interval (in wall seconds) when a restart file should be written
 * @param restart_intervals_to_keep Number of intervals that are kept in total .
 */
InputOutputManager::InputOutputManager( std::string const& input_file,
                                        std::string const& output_folder,
                                        UnitHandler const& unit_handler,
                                        OutputWriter const& output_writer,
                                        RestartManager const& restart_manager,
                                        double const time_naming_factor,
                                        std::vector<double> const& standard_output_timestamps,
                                        std::vector<double> const& interface_output_timestamps,
                                        RestoreMode const restore_mode,
                                        std::string const& restore_filename,
                                        std::vector<double> const& restart_snapshot_timestamps,
                                        int const restart_snapshot_interval,
                                        unsigned int const restart_intervals_to_keep ) :
   // Start initializer list
   unit_handler_( unit_handler ),
   logger_( LogWriter::Instance() ),
   output_writer_( output_writer ),
   restart_manager_( restart_manager ),
   output_folder_name_( output_folder ),
   time_naming_factor_( time_naming_factor ),
   standard_output_enabled_( !standard_output_timestamps.empty() ),
   standard_output_timestamps_( standard_output_timestamps ),
   interface_output_enabled_( !interface_output_timestamps.empty() ),
   interface_output_timestamps_( interface_output_timestamps ),
   restore_mode_( restore_mode ),
   restore_filename_( restore_filename ),
   restart_snapshot_timestamps_( restart_snapshot_timestamps ),
   restart_snapshot_interval_( restart_snapshot_interval ),
   restart_files_to_keep_( restart_intervals_to_keep ),
   symlink_latest_restart_name_( output_folder_name_ + RestartSubfolderName() + LatestSnapshotName() ),
   wall_time_of_last_restart_file_( std::chrono::system_clock::now() )
{
   // This Barrier is needed, otherwise we get inconsistent folder names across the ranks.
   MPI_Barrier(MPI_COMM_WORLD);
   // Only Master-rank is setting up the folders...
   if( MpiUtilities::MyRankId() == 0 ) {
      // create output folder and subfolders
      CreateOutputFolder();

      // Create the log file name
      logger_.SetLogfileName( output_folder_name_ + "/" + FileOperations::RemoveFilePath( FileOperations::RemoveFileExtension( input_file ) ) + ".log" );

      // copy inputfile into directory
      if( !input_file.empty() ) {
         std::ifstream input_stream( input_file, std::ios::binary );
         std::ofstream output_stream( output_folder_name_ + "/" + FileOperations::RemoveFilePath( input_file ), std::ios::binary );
         output_stream << input_stream.rdbuf();
         output_stream.flush();
         output_stream.close();
      }

      // Initialize the time series files (must be called after the SetupOutputFolder() to ensure the presence of folders)
      // If the time output times stamps vectors are empty no output is desired
      std::string time_series_filename;
      if( !standard_output_timestamps_.empty() ) {
         time_series_filename = OutputFileName( OutputType::Standard) + TimeSeriesSuffix();
         output_writer_.InitializeTimeSeriesFile( time_series_filename );
      }
      if( !interface_output_timestamps_.empty() ) {
         time_series_filename = OutputFileName( OutputType::Interface) + TimeSeriesSuffix();
         output_writer_.InitializeTimeSeriesFile( time_series_filename );
      }
      if constexpr( DP::DebugOutput() ) {
         time_series_filename = OutputFileName( OutputType::Debug) + TimeSeriesSuffix();
         output_writer_.InitializeTimeSeriesFile( time_series_filename );
      }
   }
}

/**
 * @brief Destructor to destroy the object with operations
 * @note In the destructor all operations are carried out to finalize the output and restart operations.
 */
InputOutputManager::~InputOutputManager() {
   // Finalizes the time series files
   std::string time_series_filename;
   if( standard_output_enabled_ ) {
      time_series_filename = OutputFileName( OutputType::Standard) + TimeSeriesSuffix();
      output_writer_.FinalizeTimeSeriesFile( time_series_filename );
   }
   if( interface_output_enabled_ ) {
      time_series_filename = OutputFileName( OutputType::Interface) + TimeSeriesSuffix();
      output_writer_.FinalizeTimeSeriesFile( time_series_filename );
   }
   if constexpr( DP::DebugOutput() ) {
      time_series_filename = OutputFileName( OutputType::Debug) + TimeSeriesSuffix();
      output_writer_.FinalizeTimeSeriesFile( time_series_filename );
   }
}

/**
 * @brief This routine creates the output folders for the simulation.
 * @note  Only the master-rank "0" is supposed to create the folder.
 */
void InputOutputManager::CreateOutputFolder() const {
   // create output folder
   FileOperations::CreateFolder( output_folder_name_ );
   // create output folder only if output is enabled
   if( standard_output_enabled_ ) {
      FileOperations::CreateFolder( output_folder_name_ + OutputSubfolderName( OutputType::Standard ) );
   }
   if( interface_output_enabled_ ) {
      FileOperations::CreateFolder( output_folder_name_ + OutputSubfolderName( OutputType::Interface ) );
   }
   if constexpr( DP::DebugOutput() ) {
      FileOperations::CreateFolder( output_folder_name_ + OutputSubfolderName( OutputType::Debug ) );
   }

   // create restart folder
   FileOperations::CreateFolder( output_folder_name_ + RestartSubfolderName() );

#ifndef PERFORMANCE
   // create an .gitignore file in the output folder
   std::ofstream output_stream( output_folder_name_ + "/.gitignore", std::ios::out );
   output_stream << "*";
   output_stream.flush();
   output_stream.close();
#endif
}

/**
 * @brief Writes all micro timestep sizes performed during the last macro timestep.
 * @param timesteps_on_finest_level Sizes of the micro timesteps.
 */
void InputOutputManager::WriteTimestepFile(  std::vector<double> const& timesteps_on_finest_level ) const {
   // can only be done for rank 0 to avoid parallel writing
   if( MpiUtilities::MyRankId() == 0 ) {
      // Create the string with all time steps
      std::string timesteps;
      for( auto const& timestep : timesteps_on_finest_level ) {
         timesteps += std::to_string( timestep ) + "\n";
      }
      // append to file
      FileOperations::AppendToTextBasedFile( output_folder_name_ + "/micro_timesteps.txt", timesteps );
   }
}

/**
 * @brief Writes the full output (all outputs desired (standard, interface, debug)) at the current timestep.
 *        If the force_output flag is set, output is written in any case.
 * @param timestep The current timestep.
 * @param force_output A flag indicating whether output should be forced.
 * @return Return whether any output was written or not.
 */
bool InputOutputManager::WriteFullOutput( double const timestep, bool const force_output ) {

   // Flag to indicate that any output has been written
   bool output_written = false;

   // dimensionalize time
   double const dimensionalized_time = unit_handler_.DimensionalizeValue( timestep, UnitType::Time );

   // string that is added to all filenames including the time
   std::string const time_name( std::to_string( dimensionalized_time * time_naming_factor_ ) );

   // Add an empty line to the logger
   logger_.LogMessage( " " );

   // Check whether if the standard output is enabled (required to prevent call .front() with empty vector )
   if( standard_output_enabled_ ) {
      // Check wether the output is forced or the next desired standard time stamp is smaller than the current
      if( force_output || standard_output_timestamps_.front() <= timestep ) {
         // erase the timestamps that are handled now (everything smaller than current time step)
         // This needs to be done since the macro timestep can jump over several given timestamps
         RemoveTimeStamps( standard_output_timestamps_, timestep );

         // Write output with logging information for standard and debug output
         std::string const output_filename( OutputFileName( OutputType::Standard) + time_name );
         std::string const time_series_filename( OutputFileName( OutputType::Standard) + TimeSeriesSuffix() );
         WriteOutput( OutputType::Standard, dimensionalized_time, output_filename, time_series_filename );

         // Set flag to true
         output_written = true;
      }
   }

   // Check whether if the interface output is enabled (required to prevent call .front() with empty vector )
   if( interface_output_enabled_ ) {
      // Check wether the output is forced or the next desired interface time stamp is smaller than the current
      if( force_output || interface_output_timestamps_.front() <= timestep ) {
         // erase the timestamps that are handled now (everything smaller than current time step)
         // This needs to be done since the macro timestep can jump over several given timestamps
         RemoveTimeStamps( interface_output_timestamps_, timestep );

         // call the output functions with the dimensionalized time
         std::string const output_filename( OutputFileName( OutputType::Interface) + time_name );
         std::string const time_series_filename( OutputFileName( OutputType::Interface) + TimeSeriesSuffix() );
         WriteOutput( OutputType::Interface, dimensionalized_time, output_filename, time_series_filename );

         // Set flag to true
         output_written = true;
      }
   }

   // Always write debug output if activated
   if constexpr( DP::DebugOutput() ) {
      std::string const output_filename( OutputFileName( OutputType::Debug) + time_name );
      std::string const time_series_filename( OutputFileName( OutputType::Debug) + TimeSeriesSuffix() );
      WriteOutput( OutputType::Debug, dimensionalized_time, output_filename, time_series_filename );

      // Set flag to true
      output_written = true;
   }

   // Add an empty line to the logger
   logger_.LogMessage( " " );
   // Flush the logger to the terminal and/or file
   logger_.Flush();

   // Returns the flag whether an output has been written
   return output_written;
}

/**
 * @brief Writes debug output with the given debug key (for detailed debug output after algorithm substep)
 * @param output_key The output key indicating at which algorithm substep the output is triggered.
 * @param output_type Output type identifier that is used (standard, interface, debug). Default: Debug
 */
void InputOutputManager::WriteSingleOutput( unsigned int const output_key, OutputType const output_type ) const {
  // correct file naming
  std::string const time_name( std::to_string( double( output_key ) ) );
  std::string const output_filename( OutputFileName( output_type ) + time_name );
  // Call the output function
  WriteOutput( output_type, double( output_key ), output_filename, "" );
}

/**
 * @brief Carries out the actual output writing (call of output writer) with additional logging
 * @param output_type Output type identifier that is used (standard, interface, debug)
 * @param output_time Time of the output
 * @param filename_without_extension Filename without extension where the output is written into
 * @param time_series_filename_without_extension Filename of the time series data that is written
 */
void InputOutputManager::WriteOutput( OutputType const output_type,
                                      double const output_time,
                                      std::string const& filename_without_extension,
                                      std::string const& time_series_filename_without_extension ) const {

   // Logging for the writing time of the output
   double const write_output_start_time = MPI_Wtime();

   // Call the output writer for writing the output
   output_writer_.WriteOutput( output_type, output_time, filename_without_extension, time_series_filename_without_extension );

   // Final logging
   logger_.LogMessage(  OutputTypeToString( output_type ) + " output file written at t = "
                      + StringOperations::ToScientificNotationString( output_time, 9 ), true, true );

   // Debug loggin information for full writing process
   if( DP::Profile() ) {
      double const seconds_elapsed = MPI_Wtime() - write_output_start_time;
      std::vector<double> all_times;
      int const number_of_ranks = MpiUtilities::NumberOfRanks();
      all_times.resize( number_of_ranks );
      //communicate all individual times
      MPI_Allgather( &seconds_elapsed, 1, MPI_DOUBLE, all_times.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD );
      std::vector<double>::iterator max_time = std::max_element( all_times.begin(), all_times.end() );
      std::vector<double>::iterator min_time = std::min_element( all_times.begin(), all_times.end() );
      double const avg_time = std::accumulate( all_times.begin(), all_times.end(), 0.0 ) / number_of_ranks;
      std::string message = "Writing " + OutputTypeToString( output_type ) + " output: avg/max(rank)/min(rank) ";
      message.append( std::to_string( avg_time ) + "/" );
      message.append( std::to_string( *max_time ) + "(" + std::to_string( distance( all_times.begin(), max_time ) ) + ")/" );
      message.append( std::to_string( *min_time ) + "(" + std::to_string( distance( all_times.begin(), min_time ) ) + ")/" );
      message.append( ") seconds elapsed." );
      logger_.LogMessage( message );
   }
}

/**
 * @brief Writes a restart snapshot if necessary at the current timestep. Old snapshots are deleted depending on user configuration.
 * @param force_output A flag indicating whether restart should be forced.
 * @param timestep The current timestep.
 * @param force_output A flag indicating whether output should be forced.
 */
void InputOutputManager::WriteRestartFile( double const timestep, bool const force_output) {
   // Check restart trigger on wall clock interval
   bool snapshot_interval_triggered = false;
   // only consider interval-based snapshots if the interval is greater zero
   if( restart_snapshot_interval_ > 0 ) {
      // only carry out for rank zero
      if( MpiUtilities::MyRankId() == 0 ) {
         std::chrono::time_point<std::chrono::system_clock> const current_wall_time = std::chrono::system_clock::now();
         int const wall_seconds_since_snapshot = std::chrono::duration_cast<std::chrono::seconds>( current_wall_time - wall_time_of_last_restart_file_ ).count();
         if( restart_snapshot_interval_ <= wall_seconds_since_snapshot ) {
            wall_time_of_last_restart_file_ = current_wall_time;
            snapshot_interval_triggered = true;
         }
      }
      // distribute restart decision among all ranks
      MPI_Bcast( &snapshot_interval_triggered, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD );
   }

   // Check restart trigger on snapshot time stamps (first check on empty vector and than compare to time step. Required to prevent call of .front() on empty vector)
   bool snapshot_timestamp_triggered = !restart_snapshot_timestamps_.empty();
   if( snapshot_timestamp_triggered ) snapshot_timestamp_triggered = restart_snapshot_timestamps_.front() <= timestep;

   // carry out restart writing if the restart is triggered by forcing, flagging or timing
   if( snapshot_interval_triggered || snapshot_timestamp_triggered || force_output ) {

      // erase the timestamps that are handled now (smaller than given timestep)
      RemoveTimeStamps( restart_snapshot_timestamps_, timestep );

      // Generate the appropriate file
      std::string const time_name( std::to_string( unit_handler_.DimensionalizeValue( timestep, UnitType::Time ) * time_naming_factor_ ) );
      std::string const filename_without_extension( RestartFileName() + time_name );

      // write the actual restart file and store its name
      std::string snapshot_filename = restart_manager_.WriteRestartFile( timestep, filename_without_extension );

      // handle filesystem access only on rank zero
      if( MpiUtilities::MyRankId() == 0 ) {

         // update symbolic link to latest snapshot file
         std::remove( symlink_latest_restart_name_.c_str() );
         [[maybe_unused]] int const result_io = symlink( FileOperations::RemoveFilePath( snapshot_filename ).c_str(), symlink_latest_restart_name_.c_str() );

         // only consider non-timestamp snapshots for deletion
         if( !snapshot_timestamp_triggered ) {
            if( restart_files_written_.size() == restart_files_to_keep_ ) {
               // remove the oldest file
               std::remove( restart_files_written_.front().c_str() );
               restart_files_written_.erase( restart_files_written_.begin() );
            }
            // add the newest file to the list
            restart_files_written_.push_back( snapshot_filename );
         }
      }
   }
}

/**
 * @brief Restores the simulation from the user configured snapshot and returns the simulation time of the snapshot.
 * @return The time of the simulation snapshot (negative if no restart is carried out)
 */
double InputOutputManager::RestoreSimulationFromSnapshot() {
   // Restart time that is returned (negative if restart is not carried out)
   double restart_time;

   // Add an empty line to the logger
   logger_.LogMessage( " " );
   // Change behavior dependent on given restore mode
   switch( restore_mode_ ) {
      case RestoreMode::Off : {
         logger_.LogMessage( "Restore from snapshot disabled!" );
         restart_time = -1.0;
      }
      break;
      case RestoreMode::Soft : {
         logger_.LogMessage( "Soft restore enabled: " );
         if( CheckIfRestoreFileExists() ) {
            logger_.LogMessage( "Initializing simulation from restart file: " + restore_filename_ );
            restart_time = restart_manager_.RestoreSimulation( restore_filename_ );
         } else {
            logger_.LogMessage( "WARNING: restart file specified in the input file does not exist!" );
            restart_time = -1.0;
         }
      }
      break;
      case RestoreMode::Forced : {
         if( !CheckIfRestoreFileExists() ) {
            throw std::runtime_error( "Force restore not possible: restart file specified in the input file does not exist" );
         }

         logger_.LogMessage( "Initializing simulation from restart file: " + restore_filename_ );
         restart_time = restart_manager_.RestoreSimulation( restore_filename_ );
      }
      break;
      default :
         // as this method is only called once at startup, we don't need to make use of the PERFORMANCE flag
         throw std::invalid_argument( "This restore mode is not known!" );
   }
   // Add an empty line to the logger
   logger_.LogMessage( " " );

   // adapt timestamp lists to new start time
   // Delete all timestamps smaller than restart time (if it is negative it will have no affect since negative time stamps are already eliminated)
   RemoveTimeStamps( standard_output_timestamps_, restart_time );
   RemoveTimeStamps( interface_output_timestamps_, restart_time );
   RemoveTimeStamps( restart_snapshot_timestamps_, restart_time );

   return restart_time;
}