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
#ifndef INSTANTIATION_INITIAL_CONDITION_H
#define INSTANTIATION_INITIAL_CONDITION_H

#include "topology/node_id_type.h"
#include "topology/id_information.h"
#include "topology/topology_manager.h"
#include "topology/tree.h"
#include "materials/material_manager.h"

#include "input_output/input_reader.h"
#include "input_output/input_reader/initial_condition_reader/initial_condition_reader.h"
#include "initial_condition/levelset_initializer.h"
#include "initial_condition/initial_condition.h"

/**
 * @brief Defines all instantiation functions required for the initial condition.
 */
namespace Instantiation {

   // factory functions
   std::vector<std::string> GetMaterialInitialConditions( InitialConditionReader const& initial_condition_reader, unsigned int number_of_materials );

   std::array<ParametricVariable, 2> CreateParametricVariables( InitialConditionReader const& initial_condition_reader );

   std::unique_ptr<LevelsetInitializer> InstantiateLevelsetInitializer( InitialConditionReader const& initial_condition_reader,
                                                                        unsigned int const levelset_index,
                                                                        std::vector<MaterialName> const& material_names,
                                                                        unsigned int const number_of_materials,
                                                                        double const node_size_on_level_zero_,
                                                                        unsigned int const maximum_level );

   // Initialization function for the initial condition class
   std::unique_ptr<InitialCondition> InstantiateInitialCondition( InputReader const& input_reader,
                                                                  TopologyManager const& topology_manager,
                                                                  Tree const& tree,
                                                                  MaterialManager const& material_manager,
                                                                  UnitHandler const& unit_handler );
}// namespace Instantiation

#endif// INSTANTIATION_INITIAL_CONDITION_H
