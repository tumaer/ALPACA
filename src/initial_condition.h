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
#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

#include <vector>
#include <array>
#include <string>
#include <memory>

#include "unit_handler.h"
#include "user_specifications/compile_time_constants.h"
#include "block_definitions/field_material_definitions.h"
#include "materials/material_definitions.h"

#include "user_expression.h"

/**
 * @brief The InitialCondition class is used to set the state of all cells according to the user input at the beginning of the simulation.
 * @note Uses the C++ Mathematical Expression Toolkit Library by Arash Partow, see respective files for License and Copyright information.
 */
class InitialCondition {
   // Instance for dimensionalization and non-dimensionalization
   UnitHandler const& unit_handler_;

   // Member variables providing the user expression of input data
   std::vector<std::string> const material_initial_expressions_;
   std::vector<std::string> const levelset_initial_expressions_;

   // Variables to define the name of parameters in the input expression
   std::vector<MaterialName> const material_names_;
   std::vector<std::string> const variable_names_prime_states_;
   std::string const variable_name_levelset_;
   std::string const variable_name_x_ = "x";
   std::string const variable_name_y_ = "y";
   std::string const variable_name_z_ = "z";
   // Additional required variables
   double const dimensionalized_node_size_on_level_zero_;
   unsigned int const maximum_level_;


   // local function to create expression that can be evaluated
   std::unique_ptr<UserExpression const> CreateInputExpression( std::string const& expression, std::vector<std::string> const& variables_out, double &x, double &y, double &z ) const;

public:
   InitialCondition() = delete;
   explicit InitialCondition( std::vector<std::string> const& material_initial_expressions,
                              std::vector<std::string> const& levelset_initial_expressions,
                              std::vector<MaterialName> const& material_names,
                              std::vector<std::string> const& variable_names_prime_states,
                              std::string const& variable_name_levelset,
                              double const node_size_on_level_zero_,
                              unsigned int const maximum_level,
                              UnitHandler const& unit_handler );
   ~InitialCondition() = default;
   InitialCondition( InitialCondition const& ) = delete;
   InitialCondition& operator=( InitialCondition const& ) = delete;
   InitialCondition( InitialCondition&& ) = delete;
   InitialCondition& operator=( InitialCondition&& ) = delete;

   // Fills the prime state buffer with appropriate values of the input expression
   void GetInitialPrimeStates( std::uint64_t const node_id, MaterialName const material, double (&initial_values)[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()] ) const;
   // Fills the levelset buffer with appropriate value of the input expressions
   void GetInitialLevelset( std::uint64_t const node_id, double (&initial_levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] ) const;
   // Gives the initial materials present
   std::vector<MaterialName> GetInitialMaterials( std::uint64_t const node_id ) const;
};

#endif // INITIAL_CONDITION_H
