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
#ifndef BOUNDARY_CONDITION_READER_H
#define BOUNDARY_CONDITION_READER_H

#include <array>

#include "enums/direction_definition.h"
#include "block_definitions/field_material_definitions.h"
#include "boundary_condition/boundary_specifications.h"

/**
 * @brief Defines the class that provides access to the boundary condition data in the input file.
 *        It serves as a proxy class for different boundary condition reader types (xml,...) that only read the actual data.
 *        Here, consistency checks are done that all read data are valid.
 */
class BoundaryConditionReader {

protected:
   // constructor can only be called from derived classes
   explicit BoundaryConditionReader() = default;

   // Functions that must be implemented by the derived classes
   virtual std::string DoReadMaterialBoundaryType( BoundaryLocation const location ) const                                        = 0;
   virtual std::string DoReadLevelSetBoundaryType( BoundaryLocation const location ) const                                        = 0;
   virtual double DoReadMaterialFixedValueBoundaryCondition( BoundaryLocation const location, std::string const& variable ) const = 0;

public:
   virtual ~BoundaryConditionReader()                        = default;
   BoundaryConditionReader( BoundaryConditionReader const& ) = delete;
   BoundaryConditionReader& operator=( BoundaryConditionReader const& ) = delete;
   BoundaryConditionReader( BoundaryConditionReader&& )                 = delete;
   BoundaryConditionReader& operator=( BoundaryConditionReader&& ) = delete;

   // Functions that must be implemented by the derived classes
   TEST_VIRTUAL MaterialBoundaryType ReadMaterialBoundaryType( BoundaryLocation const location ) const;
   TEST_VIRTUAL LevelSetBoundaryType ReadLevelsetBoundaryType( BoundaryLocation const location ) const;
   TEST_VIRTUAL std::array<double, MF::ANOP()> ReadMaterialFixedValueBoundaryConditions( BoundaryLocation const location ) const;
};

#endif// BOUNDARY_CONDITION_READER_H
