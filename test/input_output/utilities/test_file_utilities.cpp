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
#include <catch2/catch.hpp>

#include "input_output/utilities/file_utilities.h"

SCENARIO( "File extensions can be changed", "[1rank]" ) {
   GIVEN( "A filename with absolute path and .txt extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.txt";
      REQUIRE( FileUtilities::ChangeFileExtension( txt_filename_with_path, ".oop" ) == "/scratch/directory/file.oop" );
   }
   GIVEN( "A filename with absolute paths containing dots" ) {
      std::string const txt_filename_with_path = "/scratch/directory/fi.le.txt";
      REQUIRE( FileUtilities::ChangeFileExtension( txt_filename_with_path, ".oop" ) == "/scratch/directory/fi.le.oop" );
   }
   GIVEN( "A filename with long extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.longext";
      REQUIRE( FileUtilities::ChangeFileExtension( txt_filename_with_path, ".oop" ) == "/scratch/directory/file.oop" );
   }
}

SCENARIO( "File extensions can be removed", "[1rank]" ) {
   GIVEN( "A filename with absolute path and .txt extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.txt";
      REQUIRE( FileUtilities::RemoveFileExtension( txt_filename_with_path ) == "/scratch/directory/file" );
   }
   GIVEN( "A filename with absolute paths containing dots" ) {
      std::string const txt_filename_with_path = "/scratch/directory/fi.le.txt";
      REQUIRE( FileUtilities::RemoveFileExtension( txt_filename_with_path ) == "/scratch/directory/fi.le" );
   }
   GIVEN( "A filename with long extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.longext";
      REQUIRE( FileUtilities::RemoveFileExtension( txt_filename_with_path ) == "/scratch/directory/file" );
   }
}

SCENARIO( "File extensions can be obtained", "[1rank]" ) {
   GIVEN( "A filename with absolute path and .txt extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.txt";
      REQUIRE( FileUtilities::GetFileExtension( txt_filename_with_path ) == "txt" );
   }
   GIVEN( "A filename with absolute paths containing dots" ) {
      std::string const txt_filename_with_path = "/scratch/directory/fi.le.txt";
      REQUIRE( FileUtilities::GetFileExtension( txt_filename_with_path ) == "txt" );
   }
   GIVEN( "A filename with long extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.longext";
      REQUIRE( FileUtilities::GetFileExtension( txt_filename_with_path ) == "longext" );
   }
}

SCENARIO( "File paths can be removed", "[1rank]" ) {
   GIVEN( "A filename with absolute path and .txt extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.txt";
      REQUIRE( FileUtilities::RemoveFilePath( txt_filename_with_path ) == "file.txt" );
   }
   GIVEN( "A filename with absolute paths containing dots" ) {
      std::string const txt_filename_with_path = "/scratch/directory/fi.le.txt";
      REQUIRE( FileUtilities::RemoveFilePath( txt_filename_with_path ) == "fi.le.txt" );
   }
   GIVEN( "A filename with long extension" ) {
      std::string const txt_filename_with_path = "/scratch/directory/file.longext";
      REQUIRE( FileUtilities::RemoveFilePath( txt_filename_with_path ) == "file.longext" );
   }
}
