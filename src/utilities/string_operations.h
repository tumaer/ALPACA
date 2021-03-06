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
#ifndef STRING_OPERATIONS_H
#define STRING_OPERATIONS_H

#include <string>
#include <limits>

#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <type_traits>

/**
 * @brief Provide the functionality to manipulate given strings
 */
namespace StringOperations {

   // converts a number to scientific notation string, e.g. -1.0E+10
   std::string ToScientificNotationString( double const number, int const precision = std::numeric_limits<double>::digits10 + 1, bool consider_sign = false );
   // converts a string to a string without spaces and all letters are upper cases
   std::string ToUpperCaseWithoutSpaces( std::string const& word );
   // removes all spaces from a string
   std::string RemoveSpaces( std::string const& word );
   // gives a empty string with a certain width
   std::string Indent( unsigned int const width );
   // remove leading numbers from a string
   std::string RemoveLeadingNumbers( std::string const& word );
   // trim a string
   std::string Trim( std::string const& word );

   /**
    * @brief Converts a string into a variable of type T.
    * @param input The string that should be converted.
    * @tparam T Type of value to be returned.
    * @return The converted value.
    */
   template<typename T>
   T ConvertStringToValue( std::string const input ) {
      T return_value;
      std::stringstream temp( input );
      temp >> return_value;
      return return_value;
   }

   /**
    * @brief Converts string with white-space delimiter into vector of type T.
    * @param input_data String to be converted to vector.
    * @tparam T Type of vector to be returned.
    * @return Vector of values.
    */
   template<typename T>
   std::vector<T> ConvertStringToVector( std::string const input_data ) {

      //cut string at white spaces
      std::istringstream string_iterator( input_data );
      std::vector<std::string> const string_vector( std::istream_iterator<std::string>{ string_iterator },
                                                    std::istream_iterator<std::string>() );

      //return vector if type is string, otherwise convert to requested type and then return
      if constexpr( std::is_same<T, std::string>::value ) {
         return string_vector;
      } else {
         std::vector<T> return_vector( string_vector.size() );
         std::transform( std::cbegin( string_vector ), std::cend( string_vector ), std::begin( return_vector ),
                         []( std::string const& str ) { return ConvertStringToValue<T>( str ); } );
         return return_vector;
      }
   }
}// namespace StringOperations

#endif// STRING_OPERATIONS_H
