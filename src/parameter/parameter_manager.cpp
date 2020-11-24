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
#include "parameter/parameter_manager.h"

#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief Standard constructor to create the ParameterManager object.
 * @param material_manager Instance for handling models of material proeprties (such as viscosity).
 * @param halo_manager Instance to provide halo updates.
 */
ParameterManager::ParameterManager( MaterialManager const& material_manager, HaloManager& halo_manager ) : material_manager_( material_manager ),
                                                                                                           ghost_fluid_extender_( material_manager, halo_manager ) {
   /** Empty besides initializer list */
}

/**
 * @brief Updates all parameters of a single node.
 * @param node The node under consideration (indirect return).
 */
void ParameterManager::UpdateParameters( Node& node ) const {

   // Obtain the cell size for the given node
   double const cell_size = node.GetCellSize();

   // Update all parameters acting on a single material
   for( auto& phase : node.GetPhases() ) {

      // Obtain the sign of the material and the material
      std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial( phase.first );
      Material const& material        = material_manager_.GetMaterial( phase.first );

      // Start block to update the material property models
      // Shear viscosity
      if constexpr( CC::ViscosityIsActive() && CC::ShearViscosityModelActive() ) {
         // Call appropriate function depending on presence of an interface of this node
         if( node.HasLevelset() ) {
            material.GetShearViscosityModel().UpdateParameter( phase.second, cell_size, node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>(), material_sign );
         } else {
            material.GetShearViscosityModel().UpdateParameter( phase.second, cell_size );
         }
      }

      // Thermal conductivity
      if constexpr( CC::HeatConductionActive() && CC::ThermalConductivityModelActive() ) {
         // Call appropriate function depending on presence of an interface of this node
         if( node.HasLevelset() ) {
            material.GetThermalConductivityModel().UpdateParameter( phase.second, cell_size, node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>(), material_sign );
         } else {
            material.GetThermalConductivityModel().UpdateParameter( phase.second, cell_size );
         }
      }
   }
}

/**
 * @brief Extends for all nodes the parameters into the extension band.
 * @param nodes Nodes on which the extension should be performed.
 */
void ParameterManager::ExtendParameters( std::vector<std::reference_wrapper<Node>> const& nodes ) const {
   ghost_fluid_extender_.Extend( nodes );
}
