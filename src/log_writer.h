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
#ifndef LOG_WRITER_H
#define LOG_WRITER_H

#include <string>
#include <vector>

/**
 * @brief A light-weight logger to write output to the terminal and to a log file. Supports MPI execution; can log all or just the root rank.
 * @note Singleton.
 */
class LogWriter {

   std::string log_;
   std::string logfile_name_;
   int const rank_;
   bool const save_all_ranks_;

   std::string delayed_log_;

   //Singleton has only privat Constructor
   explicit LogWriter( bool const save_all_ranks );
   std::vector<std::string> FormatMessage( std::string const& message ) const;

   int ForwardRankId() const;

public:
   //Singelton "Constructor":
   static LogWriter& Instance( bool const save_all_ranks = false );

   //Singeltons may never call these methods.
   LogWriter()                   = delete;
   ~LogWriter()                  = default;
   LogWriter( LogWriter const& ) = delete;
   LogWriter& operator=( LogWriter const& ) = delete;
   LogWriter( LogWriter&& )                 = delete;
   LogWriter& operator=( LogWriter&& ) = delete;

   void SetLogfileName( std::string const name );
   void LogMessage( std::string const& message, bool const print_to_terminal = true, bool const save_in_logfile = true );
   void LogLinebreakMessage( std::string const& message, bool const print_to_terminal = true, bool const save_in_logfile = true );
   void AddBreakLine( bool const print_to_terminal = false );
   void FlushWelcomeMessage();
   void Flush();
   void FlushAlpaca( double const percentage, bool const fast_forward = false );
   void AppendDelayedLog( std::string const delayed_log );
   void DelayedLogMessage( bool const print_to_terminal = true, bool const save_in_logfile = true );
};

#endif// LOG_WRITER_H
