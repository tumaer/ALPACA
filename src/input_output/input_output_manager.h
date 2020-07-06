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
#ifndef INPUT_OUTPUT_MANAGER_H
#define INPUT_OUTPUT_MANAGER_H

#include <chrono>
#include <memory>

#include "log_writer.h"
// #include "topology/topology_manager.h"
// #include "topology/tree.h"
#include "utilities/file_operations.h"
#include "unit_handler.h"
// #include "materials/material_manager.h"
#include "input_output/output_writer.h"
#include "input_output/restart_manager.h"
#include "input_output/input_reader/input_definitions.h"
#include "input_output/output_writer/output_definitions.h"
#include "input_output/restart_manager/restart_definitions.h"

/**
 * @brief The InputOutputManager class handles creation of and access to a unique output folder
 *        and delegates all output calls. It decides whether simulation output or restart snapshots
 *        have to be written based on user configuration and calls the respective routines. Furthermore,
 *        all used micro time steps used in the simulation can be written to a file.
 */
class InputOutputManager {
   // Unit handler for dimensionalization of time
   // (keep it at this position to ensure that it it set before other are initialized that use it)
   UnitHandler const& unit_handler_;
   // Logger must not be const (otherwise no logbook cannot be appended
   LogWriter& logger_;
   // Writer for output data
   OutputWriter const& output_writer_;
   // Writer of restart data
   RestartManager const& restart_manager_;

   // Path data for output (must be first defined for initializer list in constructor)
   std::string const output_folder_name_;

   // naming factor for the time in the output files (e.g., factor = 2e4 => t=1e-8 is turned into 2e-4)
   double const time_naming_factor_;

   // Vector storing all output time stamps used in the simulation (cannot be const since already used time stamps are removed)
   bool const standard_output_enabled_;
   std::vector<double> standard_output_timestamps_;
   bool const interface_output_enabled_;
   std::vector<double> interface_output_timestamps_;

   // Restart output (vectors cannot be const due to subsequent removing operations)
   RestoreMode const restore_mode_;
   std::string const restore_filename_;
   std::vector<double> restart_snapshot_timestamps_;
   int const restart_snapshot_interval_;
   std::vector<std::string> restart_files_written_;
   unsigned int const restart_files_to_keep_;
   std::string const symlink_latest_restart_name_;
   std::chrono::time_point<std::chrono::system_clock> wall_time_of_last_restart_file_;

   // Local function to create the  appropriate folder structure
   void CreateOutputFolder() const;
   // local function to write an output with additional logging
   void WriteOutput( OutputType const output_type, double const output_time, std::string const& filename_without_extension,
                     std::string const& time_series_filename_without_extension ) const;

   // Functions for naming of files and folders
   /**
    * @brief Returns the output subfolder name.
    * @param output_type Output type to be used (standard, interface, debug).
    * @return subfolder name.
    */
   inline std::string OutputSubfolderName( OutputType const output_type ) const {
      switch( output_type ) {
         case OutputType::Debug     : { return "/debug"; }
         case OutputType::Interface : { return "/interface"; }
         default                    : { return "/domain"; } // standard output
      }
   }

   /**
    * @brief Returns the file name for the output (without time or time_series appendix).
    * @param output_type Output type to be used (standard, interface, debug).
    * @return standard output file name.
    */
   inline std::string OutputFileName( OutputType const output_type ) const {
      switch( output_type ) {
         case OutputType::Debug     : { return output_folder_name_ + OutputSubfolderName( OutputType::Debug) + "/debug_"; }
         case OutputType::Interface : { return output_folder_name_ + OutputSubfolderName( OutputType::Interface) + "/interface_"; }
         default                    : { return output_folder_name_ + OutputSubfolderName( OutputType::Standard) + "/data_"; } // standard output
      }
   }

   /**
    * @brief Returns the suffix for the time series file.
    * @return suffix of the time series file.
    */
   inline std::string TimeSeriesSuffix() const { return "time_series"; }

   /**
    * @brief Returns the restart subfolder name.
    * @return subfolder name for the restart files.
    */
   inline std::string RestartSubfolderName() const { return "/restart"; }

   /**
    * @brief Returns the file name for the restart (without time appendix).
    * @return restart file name.
    */
   inline std::string RestartFileName() const {
      return output_folder_name_ + RestartSubfolderName() + "/restart_";
   }

   /**
    * @brief Returns the name of the latest snapshot to be used for the restart.
    * @return name of latest snapshot.
    */
   inline std::string LatestSnapshotName() const { return "/latest_restart_snapshot"; }

   /**
    * @brief Checks whether the restore file specified in the input file exists.
    * @return True is the file exists, false otherwise.
    */
   inline bool CheckIfRestoreFileExists() const { return FileOperations::CheckIfPathExists( restore_filename_ ); }

public:
   explicit InputOutputManager( std::string const& input_file,
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
                                unsigned int const restart_intervals_to_keep );
   InputOutputManager() = delete;
   ~InputOutputManager();
   InputOutputManager( InputOutputManager const& ) = delete;
   InputOutputManager& operator=( InputOutputManager const& ) = delete;
   InputOutputManager( InputOutputManager&& ) = delete;
   InputOutputManager& operator=( InputOutputManager&& ) = delete;

   /**
   * @brief Checks whether the file "ABORTFILE" exists in the output folder indicating that the simulation should be aborted.
   * @return True if "ABORTFILE" exists, false otherwise.
   */
   inline bool CheckIfAbortfileExists() const { return FileOperations::CheckIfPathExists( output_folder_name_ + "/ABORTFILE" ); }

   // Function to write the time information to a file
   void WriteTimestepFile( std::vector<double> const& timesteps_on_finest_level ) const;
   // Functions to write simulation data output
   bool WriteFullOutput( double const timestep, bool const force_output = false );
   void WriteSingleOutput( unsigned int const output_key, OutputType const output_type = OutputType::Debug ) const;
   // Functions for restart
   void WriteRestartFile( double const timestep, bool const force_output = false );
   double RestoreSimulationFromSnapshot();
};

#endif // INPUT_OUTPUT_MANAGER_H
