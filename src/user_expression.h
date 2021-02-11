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
#ifndef USER_EXPRESSION_H
#define USER_EXPRESSION_H

#include <string>
#include <vector>
#include <random>

#include <exprtk.hpp>
#include "utilities/random_number_generator.h"

template<typename T>
struct random_number_expression : public exprtk::ifunction<T> {

   RandomNumberGenerator& generator_;

public:
   random_number_expression() : exprtk::ifunction<T>( 0 ), generator_( RandomNumberGenerator::Instance() ) {}
   virtual ~random_number_expression() {}

   inline T operator()() {
      return generator_.GiveRandomNumber();
   }
};

/**
 * @brief The UserExpression class represents a user defined mathematical expression with several input and output variables.
 *        It wraps the basic functionality of the powerful expression toolkit (http://www.partow.net/programming/exprtk/index.html).
 */
class UserExpression {
   random_number_expression<double> random_number_expression_;
   exprtk::symbol_table<double> symbol_table_;
   exprtk::expression<double> expression_;

public:
   UserExpression() = delete;
   explicit UserExpression( std::string const& expression,
                            std::vector<std::string> const& variables_out,
                            std::vector<std::string> const& function_variables_names,
                            std::vector<double>& function_variables_values );
   ~UserExpression()                       = default;
   UserExpression( UserExpression const& ) = delete;
   UserExpression& operator=( UserExpression const& ) = delete;
   UserExpression( UserExpression&& )                 = delete;
   UserExpression& operator=( UserExpression&& ) = delete;

   double GetValue( std::string const variable ) const;
};

#endif// USER_EXPRESSION_H
