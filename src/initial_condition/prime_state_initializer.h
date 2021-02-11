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
#ifndef PRIME_STATE_INITIALIZER_H
#define PRIME_STATE_INITIALIZER_H

#include <vector>
#include <string>

#include "unit_handler.h"
#include "topology/node_id_type.h"
#include "materials/material_definitions.h"
#include "user_specifications/compile_time_constants.h"
#include "block_definitions/field_material_definitions.h"

/**
 * @brief The PrimeStateInitializer class allows for a prime state initialization.
 * @note Uses the C++ Mathematical Expression Toolkit Library by Arash Partow, see respective files for License and Copyright information.
 */
class PrimeStateInitializer {
   // Instance for dimensionalization and non-dimensionalization
   UnitHandler const& unit_handler_;

   // Member variable providing the user expression of input data for materials
   std::vector<std::string> const prime_state_expression_strings_;
   std::vector<std::string> const prime_state_variable_names_;
   std::vector<std::string> const spatial_variable_names_ = { "x", "y", "z" };

   // Additional required variables
   double const dimensionalized_node_size_on_level_zero_;

public:
   PrimeStateInitializer( PrimeStateInitializer const& ) = delete;
   explicit PrimeStateInitializer( std::vector<std::string> const& prime_state_expression_strings,
                                   std::vector<std::string> const& prime_state_variable_names,
                                   double const dimensionalized_node_size_on_level_zero,
                                   UnitHandler const& unit_handler );
   ~PrimeStateInitializer()       = default;
   PrimeStateInitializer& operator=( PrimeStateInitializer const& ) = delete;
   PrimeStateInitializer( PrimeStateInitializer&& )                 = delete;
   PrimeStateInitializer& operator=( PrimeStateInitializer&& ) = delete;

   // Public function that can be called from outside
   void GetInitialPrimeStates( nid_t const node_id, MaterialName const material, double ( &initial_values )[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()] ) const;
};

#endif//PRIME_STATE_INITIALIZER_H
