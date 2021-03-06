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
#ifndef RESTART_DEFINITIONS_H
#define RESTART_DEFINITIONS_H

#include <stdexcept>

#include "utilities/string_operations.h"

/**
 * @brief Identifier for the restore mode used at startup.
 */
enum class RestoreMode { Off,
                         Soft,
                         Forced };

/**
 * @brief The OutputTimesType enum defines the type of the output times to be used for the output.
 *        (Off: No output is written)
 *        (Interval: A file is written for the given intervals)
 *        (Stamps: A file is written for time stamps)
 *        (IntervalStamps: A file is written for interval and time stamps)
 */
enum class SnapshotTimesType { Off,
                               Interval,
                               Stamps,
                               IntervalStamps };

/**
 * @brief Gives the proper OutputWriter type for a given string.
 * @param file_type String that should be converted.
 * @return Output writer type.
 */
inline RestoreMode StringToRestoreMode( std::string const& type ) {
   // transform string to upper case without spaces
   std::string const type_upper_case( StringOperations::ToUpperCaseWithoutSpaces( type ) );
   // switch statements cannot be used with strings
   if( type_upper_case == "OFF" ) {
      return RestoreMode::Off;
   } else if( type_upper_case == "SOFT" ) {
      return RestoreMode::Soft;
   } else if( type_upper_case == "FORCED" ) {
      return RestoreMode::Forced;
   } else {
      throw std::logic_error( "Restore mode '" + type_upper_case + "' not known!" );
   }
}

/**
 * @brief Gives the proper OutputTimes type for a given string.
 * @param times_type String that should be converted.
 * @return Output times type.
 */
inline SnapshotTimesType StringToSnapshotTimesType( std::string const& times_type ) {
   // transform string to upper case without spaces
   std::string const type_upper_case( StringOperations::ToUpperCaseWithoutSpaces( times_type ) );
   // switch statements cannot be used with strings
   if( type_upper_case == "OFF" ) {
      return SnapshotTimesType::Off;
   } else if( type_upper_case == "INTERVAL" ) {
      return SnapshotTimesType::Interval;
   } else if( type_upper_case == "STAMPS" ) {
      return SnapshotTimesType::Stamps;
   } else if( type_upper_case == "INTERVALSTAMPS" ) {
      return SnapshotTimesType::IntervalStamps;
   } else {
      throw std::logic_error( "Restart snaposhot times type '" + type_upper_case + "' not known!" );
   }
}

#endif// RESTART_DEFINITIONS_H
