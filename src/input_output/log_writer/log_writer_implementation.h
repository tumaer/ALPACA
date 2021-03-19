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

#include <filesystem>
#include <memory>

/**
 * @brief This class provides the functions to actually write logs to terminal or file output. Logged messages are buffered internally until they are actively flushed.
 * @note This class must not be used in production code. The LogWriter class is provided for this purpose.
 */
class LogWriterImplementation {

   std::unique_ptr<std::stringstream> terminal_output_;
   std::unique_ptr<std::stringstream> file_output_;

   std::unique_ptr<std::filesystem::path> logfile_;

   std::string delayed_log_;

   void InsertMessageInAllStreams( std::string const& message );

public:
   explicit LogWriterImplementation( std::unique_ptr<std::stringstream>&& terminal_output = nullptr, std::unique_ptr<std::stringstream>&& file_output = nullptr );
   ~LogWriterImplementation()                                = default;
   LogWriterImplementation( LogWriterImplementation const& ) = delete;
   LogWriterImplementation& operator=( LogWriterImplementation const& ) = delete;
   LogWriterImplementation( LogWriterImplementation&& )                 = delete;
   LogWriterImplementation& operator=( LogWriterImplementation&& ) = delete;

   void WelcomeMessage();
   void Flush();
   void SetLogfile( std::filesystem::path const& logfile );
   void RunningAlpaca( double const percentage, bool const fast_forward );
   void LogMessage( std::string const& message );
   void BufferMessage( std::string const& message );
   void LogBufferedMessages();
   void LogBreakLine();

   std::unique_ptr<std::stringstream> SwapOutTerminalOutputStream( std::unique_ptr<std::stringstream>&& new_in_old_out );
   std::unique_ptr<std::stringstream> SwapOutFileOutputStream( std::unique_ptr<std::stringstream>&& new_in_old_out );
};
