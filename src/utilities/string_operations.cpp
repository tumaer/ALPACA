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
#include "utilities/string_operations.h"

#include <sstream>
#include <iomanip>
#include <algorithm>

namespace StringOperations {

   /**
    * @brief Gives a string representation of the given input in scientific notation.
    * @param number Floating-point value.
    * @param precision Allows to control the number of displayed decimals.
    * @param consider_sign Flag if an extra space should be considered for a possible negative sign of positive numbers.
    * @return String-converted value.
    */
   std::string ToScientificNotationString( double const number, int const precision, bool consider_sign ) {
      // If number is positive add extra space for negative sign
      std::string sign_string = "";
      if( consider_sign && number >= 0.0 ) {
         sign_string = " ";
      }
      std::ostringstream out;
      out << sign_string << std::scientific << std::setprecision( precision ) << number;
      return out.str();
   }

   /**
    * @brief Converts a string into an upper-case string without white spaces.
    * @param word String that has to be converted.
    * @return converted string.
    */
   std::string ToUpperCaseWithoutSpaces( std::string const& word ) {
      // Make it local due to erasing
      std::string upper_case_word( word );
      // Remove white spaces
      upper_case_word.erase( std::remove_if( upper_case_word.begin(), upper_case_word.end(), ::isspace ), upper_case_word.end() );
      // Convert to upper case
      std::transform( upper_case_word.begin(), upper_case_word.end(), upper_case_word.begin(), ::toupper );

      return upper_case_word;
   }

   /**
    * @brief Converts a string into a string without white spaces.
    * @param word String that has to be converted.
    * @return converted string.
    */
   std::string RemoveSpaces( std::string const& word ) {
      // Make it local due to erasing
      std::string word_without_spaces( word );
      // Remove white spaces
      word_without_spaces.erase( std::remove_if( word_without_spaces.begin(), word_without_spaces.end(), ::isspace ), word_without_spaces.end() );

      return word_without_spaces;
   }

   /**
    * @brief Removes all leading and trailing whitespaces, newline characters and tabs from a string.
    * @param word String that has to be converted.
    * @return trimmed string.
    */
   std::string Trim( std::string const& word ) {
      // Make it local due to erasing
      std::string word_without_spaces( word );
      // Get the first element that does not contain any of the deilimiters
      auto const substring_start = word.find_first_not_of( " \n\t" );
      if( substring_start == std::string::npos ) {
         return "";
      }
      auto const substring_end = word.find_last_not_of( " \n\t" );
      return word.substr( substring_start, substring_end - substring_start + 1 );
   }

   /**
    * @brief Gives an empty string of certain width.
    * @param width Width of the string.
    * @return String.
    */
   std::string Indent( unsigned int const width ) {
      return std::string( width, ' ' );
   }

   /**
    * @brief Removes the leading numbers of as string.
    * @param word The string that should be modified.
    * @return The modified string.
    */
   std::string RemoveLeadingNumbers( std::string const& word ) {
      std::string name_without_leading_number( word );
      return name_without_leading_number.erase( 0, std::min( word.find_first_not_of( "0123456789" ), word.size() - 1 ) );
   }
}// namespace StringOperations
