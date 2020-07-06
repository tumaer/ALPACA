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
#ifndef PARAMETER_MANAGER_H
#define PARAMETER_MANAGER_H

#include "parameter/material_parameter_model.h"
#include "materials/material_manager.h"
#include "topology/node.h"
#include "levelset/multi_phase_manager/ghost_fluid_extender/ghost_fluid_extender_setup.h"

using GhostFluidExtenderConcretizationParameter = GhostFluidExtenderSetup::Concretize<extender>::type_parameters;

/**
 * @brief The ParameterManager handles all updates regarding the parameter buffers that lie on a single material block. In general, it is not restricted to
 *        models, such as viscosity or thermal conductivity models, that act on material properties. It can handle all models that act on a specific parameter
 *        buffer.
 */
class ParameterManager {
   // Instance for handling parameter calculation of material data (e.g., viscosity)
   MaterialManager const& material_manager_;
   // Instance to extend parameters
   GhostFluidExtenderConcretizationParameter const ghost_fluid_extender_;

public:
   ParameterManager() = delete;
   explicit ParameterManager( MaterialManager const& material_manager, HaloManager & halo_manager );
   ~ParameterManager() = default;
   ParameterManager( ParameterManager const& ) = delete;
   ParameterManager& operator=( ParameterManager const& ) = delete;
   ParameterManager( ParameterManager&& ) = delete;
   ParameterManager& operator=( ParameterManager&& ) = delete;

   // Update functions for all parameters
   void UpdateParameters( Node& node ) const;
   // Extension function for all parameters
   void ExtendParameters( std::vector<std::reference_wrapper<Node>> const& nodes ) const;
};

#endif // PARAMETER_MANAGER_H