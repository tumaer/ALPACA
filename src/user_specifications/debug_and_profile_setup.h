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
#ifndef DEBUG_PROFILE_SETUP_H
#define DEBUG_PROFILE_SETUP_H

/**
 * @brief Small class to control compile-time settings related with additional logging and output for profiling or debugging reasons.
 */
class DebugProfileSetup {

   static constexpr bool debug_logging_ = false;// Writes debug information to the logger
   static constexpr bool debug_output_  = false;// Writes the debug output (mesh and simulation data) to file
   static constexpr bool profiling_     = false;

public:
   /**
    * @brief Gives a bool indicating if any debugging option is set.
    * @return True if any debugging option is enabled. False otherwise.
    */
   static constexpr bool Debug() { return debug_logging_ || debug_output_; }

   /**
    * @brief Gives a bool to decide if additional debugging information should be provided. This function is not about debugging information
    *        for e.g. GNU Debugger, but to print additional information during a debug run or to create additional information after separate functions. All information is
    *        written to the logger class.
    *        $Convenience function. The compiler should recognize unset debug if clauses and take them out in release-builds (hopefully).$
    * @return Debug logging decision for the current build.
    */
   static constexpr bool DebugLog() { return debug_logging_; }

   /**
    * @brief Gives a bool to decide if additional debug output files are to be written out. The output is the actual mesh and simulation data of all cells that are written
    *        into the output format specified through the input file (e.g., XDMF/HDF5)
    * @return Debug output decision for the current build.
    */
   static constexpr bool DebugOutput() { return debug_output_; }

   /**
    * @brief Give a bool to indicate whether or not profiling (timing) runs are beeing run.
    * @return Profiling decision.
    */
   static constexpr bool Profile() { return profiling_; }
};

using DP = DebugProfileSetup;

#endif// DEBUG_PROFILE_SETUP_H
