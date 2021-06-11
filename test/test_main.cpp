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
#define APPROVALS_CATCH_EXISTING_MAIN
#include <ApprovalTests.hpp>
#include <mpi.h>
#ifdef __APPLE__
#include <xmmintrin.h>
#endif
#include <tuple>
#include <vector>
#include <regex>
#include <algorithm>

using namespace ApprovalTests;
auto directoryDisposer = Approvals::useApprovalsSubdirectory( "approval_tests" );

std::pair<int, int> StartUpMpiAndReportRankAndSize( int argc, char* argv[] ) {
   MPI_Init( &argc, &argv );
   int rank = -1;
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   int number_of_ranks = -1;
   MPI_Comm_size( MPI_COMM_WORLD, &number_of_ranks );
   return std::make_pair( rank, number_of_ranks );
}

/**
 * @brief Prints the starting message of the unit tests.
 * @param rank Current rank.
 */
void PrintStartupMessage( int const rank ) {
   if( rank == 0 ) {
      std::cout << "\n"
                << "                                  \\\\\n"
                << "                                  l '>\n"
                << "                                  | |\n"
                << "                                  | |\n"
                << "                                  |   Paco~\n"
                << "                                  ||    || \n"
                << "                                  ''    '' \n"
                << "                                     \n"
                << "          THE AGE OF BUGS HAS GONE - THE AGE OF UNIT-TESTS HAS COME\n"
                << "\n";
   }
}

/**
 * @brief Function that checks if the given user arguments are valid for the unit tests
 * @param argc, argv User command line arguments
 * @param number_of_ranks Number of ranks that are used
 * @return Pair with validity of tags and error message if applicable
 */
std::pair<bool, std::string> ProvidedTagsAreValid( int argc, char* argv[], int const number_of_ranks ) {
   std::vector<std::string> arguments( argv + 1, argv + argc );
   arguments.erase( std::remove_if( std::begin( arguments ), std::end( arguments ),
                                    []( std::string const& s ) { return !std::regex_search( s, std::regex( "\\[.*\\]" ) ); } ),
                    std::end( arguments ) );
   if( arguments.size() == 0 ) {
      return std::make_pair( false, "Paco needs a tag argument to run the correct suite of tests" );
   }
   if( arguments.size() >= 2 ) {
      return std::make_pair( false, "Paco needs exactly one tag argument to run the correct suite of tests" );
   }

   if( ( arguments.front() == "[1rank]" or arguments.front() == "[.slow1rank]" ) && number_of_ranks != 1 ) {
      return std::make_pair( false, "Single-core suite can only be run with one MPI-rank" );
   }

   if( arguments.front() == "[2rank]" && number_of_ranks != 2 ) {
      return std::make_pair( false, "Two-rank suite can only be run with two MPI-ranks" );
   }
   return std::make_pair( true, "" );
}

/**
 * @brief Displays the final error message to the command line.
 * @param error_message Message that is printed.
 * @param rank Current rank.
 */
void EndMpiAndDisplayErrorMessage( std::string const error_message, int const rank ) {
   if( rank == 0 ) {
      std::cout << "===============================================================================\n";
      std::cout << error_message << "\n\n";
   }
   MPI_Finalize();
}

int main( int argc, char* argv[] ) {

   auto const [rank, number_of_ranks] = StartUpMpiAndReportRankAndSize( argc, argv );

   PrintStartupMessage( rank );

   auto const [tags_valid, error_message] = ProvidedTagsAreValid( argc, argv, number_of_ranks );

   if( !tags_valid ) {
      EndMpiAndDisplayErrorMessage( error_message, rank );
      return EXIT_FAILURE;
   }

   int const result = Catch::Session().run( argc, argv );

   MPI_Finalize();

   return result;
}
