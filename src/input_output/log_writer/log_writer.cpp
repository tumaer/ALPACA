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

#include "input_output/log_writer/log_writer.h"

/**
 * @brief Default constructor using the provided streams to buffer the messages and reduced amount of external terminal of logfile writes.
 * @param terminal_output stream to use for terminal logging.
 * @param file_output stream to use for logfile logging.
 * @note Takes ownership over the provided streams.
 */
LogWriter::LogWriter( std::unique_ptr<std::stringstream>&& terminal_output, std::unique_ptr<std::stringstream>&& file_output ) : implementation_( std::move( terminal_output ), std::move( file_output ) ) {
   // Empty besides initializer list.
}

/**
 * @brief This function is used to get the Logger. If no Logger exists yet it is created, otherwise the existing logger is passed back. "Singleton Constructor"
 * @param terminal_output stream to use for terminal logging.
 * @param file_output stream to use for logfile logging.
 * @return The logger instance.
 */
LogWriter& LogWriter::Instance( std::unique_ptr<std::stringstream>&& terminal_output, std::unique_ptr<std::stringstream>&& file_output ) {
   static LogWriter instance( std::move( terminal_output ), std::move( file_output ) );
   return instance;
}

/**
 * @brief Logs the welcome (Greeting) message.
 */
void LogWriter::WelcomeMessage() {
   implementation_.WelcomeMessage();
}

/**
 * @brief Triggers writing of logged messages to terminal.
 */
void LogWriter::FlushToTerminal() {
   implementation_.FlushToTerminal();
}

/**
 * @brief Triggers writing of logged messages to file.
 */
void LogWriter::FlushToFile() {
   implementation_.FlushToFile();
}

/**
 * @brief Triggers writing of logged messages to terminal and/or file.
 */
void LogWriter::Flush() {
   implementation_.Flush();
}

/**
 * @brief Redirects logs to the specified file and logs the conducted change.
 * @param logfile New file (path) used as logfile.
 */
void LogWriter::SetLogfile( std::filesystem::path const& logfile ) {
   implementation_.SetLogfile( logfile );
}

/**
 * @brief Logs the ascii Alpaca. The position of the Alpaca implies a progress bar.
 * @param percentage Indicator how much progress is to be indicated.
 * @param fast_forward Flag whether the progress 'jumps' to the given position.
 */
void LogWriter::RunningAlpaca( double const percentage, bool const fast_forward ) {
   implementation_.RunningAlpaca( percentage, fast_forward );
}

/**
 * @brief Logs (and formats) the provided message
 * @param message The (unformatted) message.
 */
void LogWriter::LogMessage( std::string const& message ) {
   implementation_.LogMessage( message );
}

/**
 * @brief Buffers the provided message without formatting it, so it can be logged (and formatted) later.
 * @param message The (unformatted) message.
 */
void LogWriter::BufferMessage( std::string const& message ) {
   implementation_.BufferMessage( message );
}

/**
 * @brief Logs the so far buffered messages. Thereby the buffered messages are concatenated and treated as one message.
 */
void LogWriter::LogBufferedMessages() {
   implementation_.LogBufferedMessages();
}

/**
 * @brief Logs a full line of starts to indicated a sections.
 */
void LogWriter::LogBreakLine() {
   implementation_.LogBreakLine();
}
