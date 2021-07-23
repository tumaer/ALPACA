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
#include "input_output/input_output_manager.h"
#include <mpi.h>
#include <filesystem>
#include <fenv.h>// Floating-Point raising exceptions.
#ifdef __APPLE__
#include <xmmintrin.h>
#endif

#include "instantiation/input_output/instantiation_input_reader.h"
#include "instantiation/input_output/instantiation_log_writer.h"
#include "communication/mpi_utilities.h"
#include "simulation_runner.h"

/**
 * @brief Starting function of ALPACA, called from the operating system.
 * @param argc Argument count set by your operating system
 * @param argv Input arguments, MPI settings, openMP settings, etc.
 * @return Zero if program finished correctly.
 * @note Please note, due to the usage of openMP #pragmas valgrind and such tools might produce inaccurate measurements.
 * This passes valgrind tests with disabled openMP without errors
 */
int main( int argc, char* argv[] ) {

   MPI_Init( &argc, &argv );
   //Triggers signals on floating point errors, i.e. prohibits quiet NaNs and alike
#ifdef __linux__
   feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
#endif
#ifdef __APPLE__
   _MM_SET_EXCEPTION_MASK( _MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID );
#endif

   //NH separate scope for MPI.
   {
      LogWriter& logger = Instantiation::InstantiateLogWriter( MpiUtilities::MasterRank() );

      // determine the name of the executable and write it to the logger
      std::string const executable_name( argv[0] );
      logger.LogMessage( "Using executable: " + executable_name );
      logger.Flush();

      // determine the name of the input file (default: inputfile.xml)
      std::filesystem::path const input_file( argc > 1 ? argv[1] : "inputfile.xml" );
      // Instance to provide interface to the input file/data
      InputReader const input_reader( Instantiation::InstantiateInputReader( input_file ) );

      Simulation::Run( input_reader );
      logger.Flush();
   }

   MPI_Finalize();

   return 0;
}
