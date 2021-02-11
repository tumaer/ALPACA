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
#ifndef INITIAL_CONDITION_READER_H
#define INITIAL_CONDITION_READER_H

#include <string>
#include <array>
#include <vector>
#include <tuple>

#include "initial_condition/levelset_initializer_definitions.h"
#include "initial_condition/parametric_variable.h"

/**
 * @brief Defines the class that provides access to the initial condition data in the input file.
 *        It serves as a proxy class for different initial condition reader types (xml,...) that only read the actual data.
 *        Here, consistency checks are done that all read data are valid.
 */
class InitialConditionReader {

protected:
   // constructor can only be called from derived classes
   explicit InitialConditionReader() = default;

   // Functions that must be implemented by the derived classes
   virtual std::string DoReadMaterialInitialConditions( unsigned int const material_index ) const                                                                      = 0;
   virtual std::string DoReadLevelsetInitializerType( unsigned int const material_index ) const                                                                        = 0;
   virtual std::string DoReadLevelsetInitializerInput( unsigned int const levelset_index ) const                                                                       = 0;
   virtual std::vector<std::tuple<std::string, double, double, std::uint64_t>> DoReadParametricLevelsetInitializerVariables( unsigned int const levelset_index ) const = 0;
   virtual std::array<double, 3> DoReadParametricLevelsetInitializerReferencePoint( unsigned int const levelset_index ) const                                          = 0;
   virtual std::vector<std::array<double, 6>> DoReadLevelsetInitializerBoundingBoxes( unsigned int const material_index ) const                                        = 0;

public:
   virtual ~InitialConditionReader()                       = default;
   InitialConditionReader( InitialConditionReader const& ) = delete;
   InitialConditionReader& operator=( InitialConditionReader const& ) = delete;
   InitialConditionReader( InitialConditionReader&& )                 = delete;
   InitialConditionReader& operator=( InitialConditionReader&& ) = delete;

   // return functions of the reader class
   TEST_VIRTUAL std::string ReadMaterialInitialConditions( unsigned int const material_index ) const;
   TEST_VIRTUAL LevelsetInitializerType ReadLevelsetInitializerType( unsigned int const levelset_index, LevelsetInitializerType const default_type ) const;
   TEST_VIRTUAL std::string ReadLevelsetInitializerInput( unsigned int const levelset_index ) const;
   std::vector<ParametricVariable> ReadParametricLevelsetInitializerVariables( unsigned int const levelset_index ) const;
   std::array<double, 3> ReadParametricLevelsetInitializerReferencePoint( unsigned int const levelset_index ) const;
   TEST_VIRTUAL std::vector<std::array<double, 6>> ReadLevelsetInitializerBoundingBoxes( unsigned int const material_index ) const;
};

#endif// INITIAL_CONDITION_READER_H
