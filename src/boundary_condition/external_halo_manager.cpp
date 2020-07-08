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
#include "boundary_condition/external_halo_manager.h"

/**
 * @brief Constructs the external halo manager with the given input data.
 * @param material_boundary_conditions Array with already initialized material boundary conditions (in total six, one for each side in 3D)
 *        (ownership transfer takes place).
 * @param levelset_boundary_conditions Array with already initialized levelset boundary conditions (in total six, one for each side in 3D)
 *        (ownership transfer takes place).
 * @note In case of 1D and 2D simulations the array still has its maximum size, but the not used conditions are nullptr.
 */
ExternalHaloManager::ExternalHaloManager( std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6> material_boundary_conditions,
                                          std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6> levelset_boundary_conditions ) : material_boundary_conditions_( std::move( material_boundary_conditions ) ),
                                                                                                                                           levelset_boundary_conditions_( std::move( levelset_boundary_conditions ) ) {
   /** Empty beside initializer list */
}

/**
 * @brief Performs all levelset halo updates in external boundaries.
 * @param node Node for which the halo update is done (indirect return).
 * @param buffer_type Identifier of the levelset buffer that should be updated.
 * @param loc Location of the boundary to update.
 */
void ExternalHaloManager::UpdateLevelsetExternal( Node& node, InterfaceBlockBufferType const buffer_type, BoundaryLocation const loc ) const {
   levelset_boundary_conditions_[LTI( loc )]->UpdateLevelsetExternal( node, buffer_type );
}

/**
 * @brief Performs all interface tag halo updates in external boundaries.
 * @param Node for which the halo update is done (indirect return).
 * @param loc Location of the boundary to update.
 */
void ExternalHaloManager::UpdateInterfaceTagExternal( Node& node, BoundaryLocation const loc ) const {
   levelset_boundary_conditions_[LTI( loc )]->UpdateInterfaceTagExternal( node );
}

/**
 * @brief Performs all interface tag halo updates in external boundaries.
 * @param Node for which the halo update is done (indirect return).
 * @param field_type Field identifier for material block buffers.
 * @param loc Location of the boundary to update.
 */
void ExternalHaloManager::UpdateMaterialExternal( Node& node, MaterialFieldType const field_type, BoundaryLocation const loc ) const {
   material_boundary_conditions_[LTI( loc )]->UpdateMaterialExternal( node, field_type );
}