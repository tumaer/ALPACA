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
#ifndef PARAMETRIC_VARIABLE_H
#define PARAMETRIC_VARIABLE_H

#include <string>
#include "utilities/string_operations.h"

struct ParametricVariable {
   std::string name     = "";
   std::uint64_t points = 1;
   double start         = 0.0;
   double end           = 0.0;
   double delta         = 0.0;

   ParametricVariable()                            = default;
   ~ParametricVariable()                           = default;
   ParametricVariable( ParametricVariable const& ) = default;
   ParametricVariable& operator=( ParametricVariable const& ) = default;
   ParametricVariable( ParametricVariable&& )                 = default;
   ParametricVariable& operator=( ParametricVariable&& ) = default;

   /**
    * @brief Constructs a parametric variable.
    * @param name The name of the variable.
    * @param start The start point of the variable.
    * @param end The end point of the variable.
    * @param points The number of points between start and end (inclusive both).
    */
   explicit ParametricVariable( std::string const var_name, double const var_start, double const var_end, std::uint64_t const var_points ) : name( var_name ),
                                                                                                                                             points( var_points > 0 ? var_points : 1 ),
                                                                                                                                             start( CreateStart( var_start, var_end ) ),
                                                                                                                                             end( var_end ),
                                                                                                                                             delta( CreateDelta( start, end, points ) ) {
      // Empty besides initializer list
   }

   /**
    * @brief Gives the start value and makes sanity check that it is larger than the end value.
    * @return The start value.
    */
   double CreateStart( double const start, double const end ) {
      if( end < start ) {
         throw std::invalid_argument( "For the parametric variable the start value must be below the end value." );
      }
      return start;
   }

   /**
    * @brief Gives the delta value for the given start, end and points.
    * @return The delta value.
    */
   double CreateDelta( double const start, double const end, std::uint64_t const points ) {
      if( points > 1 ) {
         double const delta = ( end - start ) / double( points - 1 );
         return delta > 0 ? delta : 0.0;
      } else {
         return 0.0;
      }
   }

   /**
    * @brief Logs the data for this variable.
    * @return The log string.
    */
   std::string GetLogData( unsigned int const indent ) const {
      std::string log_string = StringOperations::Indent( indent ) + "Name: " + name + "\n";
      log_string += StringOperations::Indent( indent + 2 ) + "Start : " + StringOperations::ToScientificNotationString( start, 2 ) + "\n";
      log_string += StringOperations::Indent( indent + 2 ) + "End   : " + StringOperations::ToScientificNotationString( end, 2 ) + "\n";
      log_string += StringOperations::Indent( indent + 2 ) + "Delta : " + StringOperations::ToScientificNotationString( delta, 2 ) + "\n";
      log_string += StringOperations::Indent( indent + 2 ) + "Points: " + std::to_string( points ) + "\n";
      return log_string;
   }
};

#endif// PARAMETRIC_VARIABLE_H
