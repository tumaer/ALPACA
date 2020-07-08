
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
#include "utilities/file_operations.h"

#include <sys/stat.h>
#include <sstream>
#include <fstream>

namespace FileOperations {

   /**
 * @brief Removes the file extension from the given filename.
 * @param filename.
 * @return The name of the file without extension.
 */
   std::string RemoveFileExtension( std::string const& filename ) {
      return filename.substr( 0, filename.find_last_of( "." ) );
   }

   /**
 * @brief Removes the path from the given filename.
 * @param filename.
 * @return The name of the file without path.
 */
   std::string RemoveFilePath( std::string const& filename ) {
      return filename.substr( filename.find_last_of( "/" ) + 1 );
   }

   /**
 * @brief Changes the extension of a given file. Can handle filenames with full (absolute or relative) path.
 * @param filename_with_path .
 * @param extension The future extension of the file.
 * @return New filename (with path if input had path).
 */
   std::string ChangeFileExtension( std::string const& filename_with_path, std::string const extension ) {
      auto const positon_of_last_dot = filename_with_path.find_last_of( '.' );
      return filename_with_path.substr( 0, positon_of_last_dot ) + extension;
   }

   /**
 * @brief Gives the extension of a given file. Can handle filenames with full (absolute or relative) path.
 * @param filename Name of the file.
 * @return Extension of the file (without dot).
 */
   std::string GetFileExtension( std::string const& filename ) {
      auto const positon_of_last_dot = filename.find_last_of( '.' );
      return filename.substr( positon_of_last_dot + 1 );
   }

   /**
 * @brief Checks if the given path already exists.
 * @param path The name of the path, whose existence is to be checked.
 * @return "true" if the path exists, otherwise false.
 */
   bool CheckIfPathExists( std::string const& path ) {
      struct stat info;
      return ( stat( path.c_str(), &info ) == 0 );
   }

   /**
 * @brief Creates a folder with the name "path".
 * @param path The name of the folder to be created.
 * @return "true" if the folder "path" was successfully created, otherwise false.
 */
   bool CreateFolder( std::string const& path ) {
      mode_t mode = 0755;
      int ret     = mkdir( path.c_str(), mode );
      return ( ret == 0 );
   }

   /**
 * @brief Finds and adds the next free integer-counter to the path.
 * $THIS FUNCTION IS NOT THREAD SAFE, i.e. if one thread creates a folder in between calls, the results are different$.
 * @param path The candidate for the path name.
 * @param start_number The next number to be checked (0 by default).
 * @return The path name that is not yet existing.
 */
   std::string AddUnusedNumberToPath( std::string const& path, unsigned int const start_number ) {

      std::string target;
      if( start_number == 0 ) {
         target = path;
      } else {
         target = path + "-" + std::to_string( start_number );
      }

      if( CheckIfPathExists( target ) ) {
         target = AddUnusedNumberToPath( path, start_number + 1 );
      }

      return target;
   }

   /**
 * @brief Writes a file containing human readable (ascii) text to disk.
 * @param filename The name of the file to be created.
 * @param content The text to be written into the file.
 */
   void WriteTextBasedFile( std::string const& filename, std::string const& content ) {
      std::ofstream output_stream( filename, std::ios::trunc );
      output_stream << content;
      output_stream.flush();
      output_stream.close();
   }

   /**
 * @brief Appends content to an existing file containing human readable (ascii) text to disk.
 * @param filename The name of the file to be appended.
 * @param content The text to be written into the file.
 */
   void AppendToTextBasedFile( std::string const& filename, std::string const& content ) {
      std::ofstream output_stream( filename, std::ios::app );
      output_stream << content;
      output_stream.flush();
      output_stream.close();
   }

}// namespace FileOperations