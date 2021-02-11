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
#ifndef EXTERNAL_HALO_MANAGER_H
#define EXTERNAL_HALO_MANAGER_H

#include <array>
#include <memory>
#include "boundary_condition/boundary_specifications.h"
#include "boundary_condition/material_boundary_condition.h"
#include "boundary_condition/levelset_boundary_condition.h"
#include "topology/node.h"

/**
 * @brief Container of the external boundaries conditions for materials and interfaces. Furthermore, provides functionality to update halo cells
 *        of external boundaries of a single node, based on the provided boundary conditions.
 */
class ExternalHaloManager {

private:
   // arrays with full initialized boundary condition on each location (east, west, north, south, top bottom). If not present it is a nullptr.
   std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6> const material_boundary_conditions_;
   std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> const levelset_boundary_conditions_;

public:
   ExternalHaloManager() = delete;
   explicit ExternalHaloManager( std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6> material_boundary_conditions,
                                 std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> levelset_boundary_conditions );
   ~ExternalHaloManager()                            = default;
   ExternalHaloManager( ExternalHaloManager const& ) = delete;
   ExternalHaloManager& operator=( ExternalHaloManager const& ) = delete;
   ExternalHaloManager( ExternalHaloManager&& )                 = delete;
   ExternalHaloManager& operator=( ExternalHaloManager&& ) = delete;

   void UpdateLevelsetExternal( Node& node, InterfaceBlockBufferType const buffer_type, BoundaryLocation const loc ) const;
   void UpdateInterfaceTagExternal( std::int8_t ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()], BoundaryLocation const loc ) const;
   void UpdateMaterialExternal( Node& node, MaterialFieldType const field_type, BoundaryLocation const loc ) const;
};

#endif /* EXTERNAL_HALO_MANAGER_H */
