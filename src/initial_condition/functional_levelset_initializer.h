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
#ifndef FUNCTIONAL_LEVELSET_INITIALIZER_H
#define FUNCTIONAL_LEVELSET_INITIALIZER_H

#include "user_expression.h"
#include "initial_condition/levelset_initializer.h"
#include "block_definitions/field_interface_definitions.h"

/**
 * @brief The FunctionalLevelsetInitializer class allows for a levelset initialization based on a levelset function. It contains the levelset computation,
 *        the loop structure is in the base class.
 * @note Uses the C++ Mathematical Expression Toolkit Library by Arash Partow, see respective files for License and Copyright information.
 * @note For the functional levelset initializer some things need to be considered:
 *       1. The function must be a signed distance function.
 *       2. The signed distance must be in the magnitude O(1) and must be in the range [-8, 8] in the first cells near the interface.
 */
class FunctionalLevelsetInitializer : public LevelsetInitializer {
   // Member variables for this class only
   std::string const levelset_variable_name_ = std::string( IF::InputName( InterfaceDescription::Levelset ) );
   std::string const levelset_expression_string_;
   std::vector<double> expression_point_ = { 0.0, 0.0, 0.0 };
   UserExpression const levelset_expression_;

   // Functions required from base class
   double ComputeSignedLevelsetValue( std::array<double, 3> const& point ) override;
   std::string GetTypeLogData( unsigned int const indent ) const override;

public:
   FunctionalLevelsetInitializer() = delete;
   explicit FunctionalLevelsetInitializer( std::string const& levelset_expression_string,
                                           std::vector<std::array<double, 6>> const& bounding_boxes,
                                           std::vector<MaterialName> const& material_names,
                                           double const node_size_on_level_zero,
                                           unsigned int const maximum_level );
   virtual ~FunctionalLevelsetInitializer()                              = default;
   FunctionalLevelsetInitializer( FunctionalLevelsetInitializer const& ) = delete;
   FunctionalLevelsetInitializer& operator=( FunctionalLevelsetInitializer const& ) = delete;
   FunctionalLevelsetInitializer( FunctionalLevelsetInitializer&& )                 = delete;
   FunctionalLevelsetInitializer& operator=( FunctionalLevelsetInitializer&& ) = delete;
};

#endif//FUNCTIONAL_LEVELSET_INITIALIZER_H
