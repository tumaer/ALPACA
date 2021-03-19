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

#include <catch2/catch.hpp>
#include <sstream>
#include <string>
#include "input_output/log_writer/log_writer_implementation.h"

/**
 * @brief Helper enum to reuse the same functionality when checking the logger's different streams.
 */
enum class OutputPlace { Terminal,
                         File };

/**
 * @brief Collection of strings expected as logger output in the test cases.
 */
namespace ExpectedStrings {

   /**
    * @brief Collection of message parts frequently contained in logger output.
    */
   namespace BuildingBlocks {
      // Some time in the future compiles will support constexpr  string, until then we stick with boring normal const
      std::string const star_line    = "|************************************************************************************|\n";
      std::string const empty_line   = "|*                                                                                  *|\n";
      std::string const ascii_alpaca = "|*  \\\\                                                                              *|\n"
                                       "|*  l '>                                                                            *|\n"
                                       "|*  | |                                                                             *|\n"
                                       "|*  | |                                                                             *|\n"
                                       "|*  | alpaca~                                                                       *|\n"
                                       "|*  ||    ||                                                                        *|\n"
                                       "|*  ''    ''                                                                        *|\n";
      std::string const greetings_message_terminal = "|*  THE AGE OF ALPACA HAS COME                                                      *|\n";
      std::string const greetings_message_file     = "|* TERMINAL OUTPUT HAS GONE - THE ALPACA LOGBOOK HAS COME                           *|\n";
      std::string const short_message              = "|* Hello World!                                                                     *|\n";
      std::string const long_message               = "|* Waltz, bad nymph, for quick jigs vex. Jived fox nymph grabs quick waltz. Glib jo *|\n"
                                       "|* cks quiz nymph to vex dwarf. Sphinx of black quartz, judge my vow. How vexingly  *|\n"
                                       "|* quick daft zebras jump!                                                          *|\n";
      std::string const breaklinen_message = "|* We break the first line here                                                     *|\n"
                                             "|* The second one there                                                             *|\n"
                                             "|* and the lastline ends without a break but just a period.                         *|\n";
      std::string const running_alpaca_zero = "|*       \\\\                                                                         *|\n"
                                              "|*       l '>                                                                       *|\n"
                                              "|*       | |                                                                        *|\n"
                                              "|*       | |                                                                        *|\n"
                                              "|* ~alpaca |                                                                        *|\n"
                                              "|*  ||    ||                                                                        *|\n"
                                              "|*  ''    ''                                                                        *|\n";

      std::string const running_alpaca_seventeen = "|*                  \\\\                                                              *|\n"
                                                   "|* >>>>>>>>>>>>>>>> l '>                                                            *|\n"
                                                   "|*                  | |                                                             *|\n"
                                                   "|* >>>>>>>>>>>>>>>> | |                                                             *|\n"
                                                   "|*            ~alpaca |                                                             *|\n"
                                                   "|* >>>>>>>>>>> ||    ||                                                             *|\n"
                                                   "|*             ''    ''                                                             *|\n";

      std::string const running_alpaca_sixtysix = "|*                                                     \\\\                           *|\n"
                                                  "|*                                                     l '>                         *|\n"
                                                  "|*                                                     | |                          *|\n"
                                                  "|*                                                     | |                          *|\n"
                                                  "|*                                               ~alpaca |                          *|\n"
                                                  "|*                                                ||    ||                          *|\n"
                                                  "|*                                                ''    ''                          *|\n";

      std::string const running_alpaca_onehundred = "|*                                                                             \\\\   *|\n"
                                                    "|* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> l '> *|\n"
                                                    "|*                                                                             | |  *|\n"
                                                    "|* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> | |  *|\n"
                                                    "|*                                                                       ~alpaca |  *|\n"
                                                    "|* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ||    ||  *|\n"
                                                    "|*                                                                        ''    ''  *|\n";

      std::string const delayed_message = "|* This message is brought to you by our sponsors:                                  *|\n"
                                          "|* ALPACA is free for all thanks to funding by <insert funding agencies>            *|\n"
                                          "|* First delayed log message. Contains breaklines here                              *|\n"
                                          "|* and here.                                                                        *|\n"
                                          "|* Fun, right?Second delayed log message. No breaklines.                            *|\n";

      std::string const sandwich_message         = "|* Sandwich Message                                                                 *|\n";
      std::string const inital_logfile_message   = "|* Logger: Changed log file from: Not Set to a.out                                  *|\n";
      std::string const exchange_logfile_message = "|* Logger: Changed log file from: a.out to b.log                                    *|\n";

   }// namespace BuildingBlocks

   /**
    * @brief The expected result of loggers welcome-message test.
    * @param op The logger's stream beeing tested.
    */
   std::string WelcomeMessage( OutputPlace const op ) {
      std::string const greetings = op == OutputPlace::Terminal ? BuildingBlocks::greetings_message_terminal : BuildingBlocks::greetings_message_file;
      return BuildingBlocks::star_line + BuildingBlocks::ascii_alpaca + BuildingBlocks::empty_line + greetings + BuildingBlocks::empty_line + BuildingBlocks::star_line;
   }

   /**
    * @brief The expected result of loggers running-alpacas test.
    */
   std::string RunningAlpacas() {
      return BuildingBlocks::star_line + BuildingBlocks::empty_line + BuildingBlocks::running_alpaca_zero + BuildingBlocks::empty_line + BuildingBlocks::star_line + BuildingBlocks::star_line + BuildingBlocks::empty_line + BuildingBlocks::running_alpaca_seventeen + BuildingBlocks::empty_line + BuildingBlocks::star_line + BuildingBlocks::star_line + BuildingBlocks::empty_line + BuildingBlocks::running_alpaca_sixtysix + BuildingBlocks::empty_line + BuildingBlocks::star_line + BuildingBlocks::star_line + BuildingBlocks::empty_line + BuildingBlocks::running_alpaca_onehundred + BuildingBlocks::empty_line + BuildingBlocks::star_line;
   }

   /**
    * @brief The expected result of loggers logged-messages test.
    */
   std::string LoggedMessages() {
      return BuildingBlocks::short_message + BuildingBlocks::long_message + BuildingBlocks::breaklinen_message;
   }

   /**
    * @brief The expected result of loggers delayed-logging test.
    */
   std::string DelayedMessages() {
      return BuildingBlocks::delayed_message;
   }

   /**
    * @brief The expected result of loggers sandwich-message test.
    */
   std::string SandwichMessage() {
      return BuildingBlocks::star_line + BuildingBlocks::sandwich_message + BuildingBlocks::star_line;
   }

   /**
    * @brief The expected result of loggers logfile-switch test.
    */
   std::string LogfileSwitch() {
      return BuildingBlocks::inital_logfile_message + BuildingBlocks::exchange_logfile_message;
   }
}// namespace ExpectedStrings

namespace {
   /**
    * @brief Helper function obtaining the message the logger would write to terminal or logfile.
    * @param logger The logger (implementation) instance under test.
    * @return The two strings residing in the terminal and logfile streams.
    */
   std::pair<std::string, std::string> LoggerOutput( LogWriterImplementation& logger ) {
      std::unique_ptr<std::stringstream> future_rerouted_output_terminal = std::make_unique<std::stringstream>( std::ios_base::out );
      std::unique_ptr<std::stringstream> future_rerouted_output_file     = std::make_unique<std::stringstream>( std::ios_base::out );
      auto const written_terminal                                        = logger.SwapOutTerminalOutputStream( std::move( future_rerouted_output_terminal ) );
      auto const written_file                                            = logger.SwapOutFileOutputStream( std::move( future_rerouted_output_file ) );
      return std::make_pair( written_terminal->str(), written_file->str() );
   }
}// namespace

SCENARIO( "Single-Rank Log Writing", "[1rank]" ) {
   GIVEN( "A log writer (implementation)" ) {
      std::unique_ptr<std::stringstream> rerouted_output_terminal = std::make_unique<std::stringstream>( std::ios_base::out );
      std::unique_ptr<std::stringstream> rerouted_output_file     = std::make_unique<std::stringstream>( std::ios_base::out );
      LogWriterImplementation logger                              = LogWriterImplementation( std::move( rerouted_output_terminal ), std::move( rerouted_output_file ) );
      WHEN( "We swap in and out the stream without flushing" ) {
         std::unique_ptr<std::stringstream> terminal_swap = std::make_unique<std::stringstream>( std::ios_base::out );
         std::unique_ptr<std::stringstream> file_swap     = std::make_unique<std::stringstream>( std::ios_base::out );
         std::string const inital_message_terminal        = "Terminal Swap";
         std::string const inital_message_file            = "File Swap";
         *terminal_swap << inital_message_terminal;
         *file_swap << inital_message_file;
         auto&& original_terminal                    = logger.SwapOutTerminalOutputStream( std::move( terminal_swap ) );
         auto&& original_file                        = logger.SwapOutFileOutputStream( std::move( file_swap ) );
         std::string const original_terminal_content = original_terminal->str();
         std::string const original_file_content     = original_file->str();
         terminal_swap                               = logger.SwapOutTerminalOutputStream( std::move( original_terminal ) );
         file_swap                                   = logger.SwapOutFileOutputStream( std::move( original_file ) );
         THEN( "The original stream contents are empty" ) {
            REQUIRE( original_terminal_content == "" );
            REQUIRE( original_file_content == "" );
         }
         THEN( "The swapped in and out streams still hold the inital message" ) {
            REQUIRE( terminal_swap->str() == inital_message_terminal );
            REQUIRE( file_swap->str() == inital_message_file );
         }
      }

      WHEN( "We log the welcome message" ) {
         logger.WelcomeMessage();
         THEN( "The welcome message with correct greeting is found in each stream" ) {
            auto const [terminal_output, file_output] = LoggerOutput( logger );
            REQUIRE( terminal_output == ExpectedStrings::WelcomeMessage( OutputPlace::Terminal ) );
            REQUIRE( file_output == ExpectedStrings::WelcomeMessage( OutputPlace::File ) );
         }
      }

      WHEN( "We log a short then a long message and finally a message containing linebreaks" ) {
         logger.LogMessage( "Hello World!" );
         logger.LogMessage( "Waltz, bad nymph, for quick jigs vex. Jived fox nymph grabs quick waltz. Glib jocks quiz nymph to vex dwarf. Sphinx of black quartz, judge my vow. How vexingly quick daft zebras jump!" );
         logger.LogMessage( "We break the first line here\nThe second one there\nand the lastline ends without a break but just a period." );
         THEN( "The streams hold three lines one containing the framed short message and two containing the framed broken-up long message" ) {
            auto const [terminal_output, file_output] = LoggerOutput( logger );
            REQUIRE( terminal_output == file_output );
            REQUIRE( terminal_output == ExpectedStrings::LoggedMessages() );
         }
      }

      WHEN( "We log four running Alpacas at 0, 17, 66 and 100 percent, with the second and last fast-forwarded" ) {
         logger.RunningAlpaca( 0, false );
         logger.RunningAlpaca( 0.17, true );
         logger.RunningAlpaca( 0.66, false );
         logger.RunningAlpaca( 1.0, true );
         THEN( "The the four (forwarded) alpacas are found in both streams" ) {
            auto const [terminal_output, file_output] = LoggerOutput( logger );
            REQUIRE( terminal_output == file_output );
            REQUIRE( terminal_output == ExpectedStrings::RunningAlpacas() );
         }
      }
      WHEN( "We log a normal message, then two delayed ones and a normal message again and finally log the dellayed messages" ) {
         logger.LogMessage( "This message is brought to you by our sponsors:" );
         logger.BufferMessage( "First delayed log message. Contains breaklines here\nand here.\nFun, right?" );
         logger.BufferMessage( "Second delayed log message. No breaklines." );
         logger.LogMessage( "ALPACA is free for all thanks to funding by <insert funding agencies>" );
         logger.LogBufferedMessages();
         THEN( "The messages appear first, followed by the delayed messages, which are squashed together" ) {
            auto const [terminal_output, file_output] = LoggerOutput( logger );
            REQUIRE( terminal_output == file_output );
            REQUIRE( terminal_output == ExpectedStrings::DelayedMessages() );
         }
      }
      WHEN( "We log a message in a break line sandwich" ) {
         logger.LogBreakLine();
         logger.LogMessage( "Sandwich Message" );
         logger.LogBreakLine();
         THEN( "The sandwich message apperase surrounded by stars" ) {
            auto const [terminal_output, file_output] = LoggerOutput( logger );
            REQUIRE( terminal_output == file_output );
            REQUIRE( terminal_output == ExpectedStrings::SandwichMessage() );
         }
      }
      WHEN( "We initally set the logfile to a.out then change it agin to b.log" ) {
         logger.SetLogfile( "a.out" );
         logger.SetLogfile( "b.log" );
         THEN( "We expect to the the switch in the log streams" ) {
            auto const [terminal_output, file_output] = LoggerOutput( logger );
            REQUIRE( terminal_output == file_output );
            REQUIRE( terminal_output == ExpectedStrings::LogfileSwitch() );
         }
      }
   }
}
