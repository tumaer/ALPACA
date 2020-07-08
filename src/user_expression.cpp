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
#include "user_expression.h"

#include <string>
#include <vector>

/**
 * @brief Creates a new UserExpression object from an expression and several input and output variables.
 * @param expression_string The expression to be compiled.
 * @param variables_in The input variable names and references to their values.
 * @param variables_out The out variable names.
 */
UserExpression::UserExpression( std::string const expression_string, std::vector<std::tuple<std::string, double&>> const variables_in, std::vector<std::string> const variables_out ) : random_number_expression_() {
   symbol_table_.add_function( "rand", random_number_expression_ );

   for( auto& var : variables_in ) {
      symbol_table_.add_variable( std::get<0>( var ), std::get<1>( var ) );
   }

   for( auto& var : variables_out ) {
      symbol_table_.create_variable( var );
   }

   symbol_table_.add_constants();

   expression_.register_symbol_table( symbol_table_ );

   exprtk::parser<double> parser;
   // there might be variables in the expression that are not registered (e.g. z velocity in 2D cases) and should be resolved automatically
   parser.enable_unknown_symbol_resolver();
   if( !parser.compile( expression_string, expression_ ) ) {
      throw std::logic_error( "Error in expression: " + parser.error() + " Expression: " + expression_string );
   }
}

/**
 * @brief Evaluates the expression and returns the value of the specified variable.
 * @param variable The name of the variable whose value should be returned.
 * @return The value of the specified variable.
 */
double UserExpression::GetValue( std::string const variable ) const {
   expression_.value();
   return symbol_table_.get_variable( variable )->value();
}
