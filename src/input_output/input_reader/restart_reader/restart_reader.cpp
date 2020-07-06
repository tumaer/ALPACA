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
#include "input_output/input_reader/restart_reader/restart_reader.h"

#include <algorithm>
#include <stdexcept>

/**
 * @brief Gives the checked mode used to restore a simulation.
 * @return restore mode identifier of the simulation.
 */
RestoreMode RestartReader::ReadRestoreMode() const {
   return StringToRestoreMode( DoReadRestoreMode() );
}

/**
 * @brief Gives the checked filename to be used for the output.
 * @return name of the restore file.
 * @note Filename without any white spaces.
 */
std::string RestartReader::ReadRestoreFilename() const {
   return StringOperations::RemoveSpaces( DoReadRestoreFilename() );
}

/**
 * @brief Gives the checked times type used for writing the restart snapshots.
 * @return snapshot times type identifier.
 */
SnapshotTimesType RestartReader::ReadSnapshotTimesType() const {
   return StringToSnapshotTimesType( DoReadSnapshotTimesType() );
}

/**
 * @brief Gives the checked interval to be used for writign restart files (in wall seconds).
 * @return interval to be used.
 */
unsigned int RestartReader::ReadSnapshotInterval() const {
   // read and make consistency check
   int const interval( DoReadSnapshotInterval() );
   if( interval < 0 ) {
      throw std::invalid_argument( "Snapshot interval for restart files must NOT be below zero!" );
   }

   return static_cast<unsigned int>( interval );
}

/**
 * @brief Gives the checked number of snapshots to be kept (for the interval based writing).
 * @return number of snapshots.
 */
unsigned int RestartReader::ReadSnapshotIntervalsToKeep() const {
   // read and make consistency check
   int const keep( DoReadSnapshotIntervalsToKeep() );
   if( keep < 0 ) {
      throw std::invalid_argument( "Snapshot intervals to keep for restart files must NOT be below zero!" );
   }

   return static_cast<unsigned int>( keep );
}

/**
 * @brief Gives the checked time stamps to be used to write restart snapshots
 * @return timestamps to be used
 */
std::vector<double> RestartReader::ReadSnapshotTimeStamps() const {
   // Obtain the time stamps and sort them before return
   std::vector<double> time_stamps( DoReadSnapshotTimeStamps() );
   // If negative elements are present throw error
   if( std::any_of( time_stamps.begin(), time_stamps.end(), [] ( double const timestamp ) { return timestamp < 0.0; } ) ) {
      throw std::invalid_argument( "All time stamps for the restart snapshots must be positive or zero!" );
   }
   // Sort the time stamps
   std::sort( time_stamps.begin(), time_stamps.end() );

   return time_stamps;
}