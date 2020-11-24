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
#ifndef TWO_PHASE_MANAGER_H
#define TWO_PHASE_MANAGER_H

#include "multi_phase_manager.h"
#include "buffer_handler.h"
#include "interface_tags/interface_tag_functions.h"

/**
 * @brief The TwoPhaseManager provides functionality to perform two-phase flow simulation by a single-level set method.
 */
class TwoPhaseManager : public MultiPhaseManager<TwoPhaseManager> {

   friend MultiPhaseManager;

private:
   template<InterfaceDescriptionBufferType T>
   void SetVolumeFractionBuffer( Node& node ) const;

   void MixImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const;
   void EnforceWellResolvedDistanceFunctionImplementation( std::vector<std::reference_wrapper<Node>> const& nodes, bool const is_last_stage = false ) const;
   void ExtendPrimeStatesImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const;
   void ExtendInterfaceStatesImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const;
   void UpdateIntegratedBufferImplementation( std::vector<std::reference_wrapper<Node>> const& nodes, bool const is_last_stage ) const;
   void PropagateLevelsetImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const;
   void InitializeVolumeFractionBufferImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const;
   void ObtainInterfaceStatesImplementation( std::vector<std::reference_wrapper<Node>> const& nodes, bool const reset_interface_states = false ) const;

   /**
    * @brief Sets the interface tags on the finest level. Implementation for the reinitialized buffer.
    * @param nodes_containing_level_set The nodes on the finest level, which have a level-set block.
    * @tparam IDB The interface buffer type for the interface tag update.
    */
   template<InterfaceDescriptionBufferType IDB>
   void UpdateInterfaceTagsOnFinestLevel( std::vector<std::reference_wrapper<Node>> const& nodes_containing_level_set ) const {

      for( Node& node : nodes_containing_level_set ) {
         InterfaceTagFunctions::SetInternalCutCellTagsFromLevelset( node.GetInterfaceBlock().GetInterfaceDescriptionBuffer<IDB>()[InterfaceDescription::Levelset], node.GetInterfaceTags<IDB>() );
      }

      halo_manager_.InterfaceTagHaloUpdateOnLmax<IDB>();

      for( Node& node : nodes_containing_level_set ) {
         InterfaceTagFunctions::SetTotalInterfaceTagsFromCutCells( node.GetInterfaceTags<IDB>() );
      }
      halo_manager_.InterfaceTagHaloUpdateOnLmax<IDB>();
   }

public:
   TwoPhaseManager() = delete;
   explicit TwoPhaseManager( MaterialManager const& material_manager, HaloManager& halo_manager );
   ~TwoPhaseManager()                        = default;
   TwoPhaseManager( TwoPhaseManager const& ) = delete;
   TwoPhaseManager& operator=( TwoPhaseManager const& ) = delete;
   TwoPhaseManager( TwoPhaseManager&& )                 = delete;
   TwoPhaseManager& operator=( TwoPhaseManager&& ) = delete;
};

#endif//TWO_PHASE_MANAGER_H
