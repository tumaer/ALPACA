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
#ifndef MULTI_PHASE_MANAGER_H
#define MULTI_PHASE_MANAGER_H

#include "user_specifications/numerical_setup.h"
#include "communication/communication_manager.h"
#include "materials/material_manager.h"
#include "interface_interaction/interface_state_calculator.h"

#include "levelset/geometry/geometry_calculator_setup.h"
#include "buffer_handler_setup.h"
#include "cut_cell_mixer/cut_cell_mixer_setup.h"
#include "levelset_reinitializer/levelset_reinitializer_setup.h"
#include "scale_separator/scale_separator_setup.h"
#include "ghost_fluid_extender/ghost_fluid_extender_setup.h"
#include "interface_extender/interface_extender_setup.h"

using GeometryCalculatorConcretization    = GeometryCalculatorSetup::Concretize<geometry_calculator>::type;
using BufferHandlerConcretization         = BufferHandlerSetup::Concretize<buffer_handler>::type;
using CutCellMixerConcretization          = CutCellMixerSetup::Concretize<cut_cell_mixer>::type;
using LevelsetReinitializerConcretization = LevelsetReinitializerSetup::Concretize<levelset_reinitializer>::type;
using ScaleSeparatorConcretization        = ScaleSeparatorSetup::Concretize<scale_separator>::type;
using GhostFluidExtenderConcretization    = GhostFluidExtenderSetup::Concretize<extender>::type_primestates;
using InterfaceExtenderConcretization     = InterfaceExtenderSetup::Concretize<interface_extender>::type_states;

/**
 * @brief The MultiPhaseManager provides functionality to simulate multi-phase flows. It allows to propagate a level-set field in time, to perform cut-cell mixing
 *        and to extend fluid states to ghost cells.
 * @tparam DerivedMultiPhaseManager Typename as template parameter due to CRTP.
 */
template<typename DerivedMultiPhaseManager>
class MultiPhaseManager {

   friend DerivedMultiPhaseManager;
   // general classes without concretization
   MaterialManager const& material_manager_;
   HaloManager& halo_manager_;// TODO-19 NH Think about making it const (rats tail)
   // class concretization use
   GeometryCalculatorConcretization const geometry_calculator_;
   BufferHandlerConcretization const buffer_handler_;
   InterfaceStateCalculator const interface_state_calculator_;
   CutCellMixerConcretization const cut_cell_mixer_;
   LevelsetReinitializerConcretization const levelset_reinitializer_;
   ScaleSeparatorConcretization const scale_separator_;
   GhostFluidExtenderConcretization const ghost_fluid_extender_;
   InterfaceExtenderConcretization const interface_extender_;

   /**
    * @brief Default constructor for a MultiPhaseManager.
    * @param material_manager Instance of a material manager, which already has been initialized according to the user input.
    * @param halo_manager Instance to a HaloManager which provides MPI-related methods.
    */
   explicit MultiPhaseManager( MaterialManager const& material_manager, HaloManager& halo_manager ) : material_manager_( material_manager ),
                                                                                                      halo_manager_( halo_manager ),
                                                                                                      geometry_calculator_(),
                                                                                                      buffer_handler_( material_manager_ ),
                                                                                                      interface_state_calculator_( material_manager_ ),
                                                                                                      cut_cell_mixer_( halo_manager, material_manager_ ),
                                                                                                      levelset_reinitializer_( halo_manager ),
                                                                                                      scale_separator_( material_manager_, halo_manager ),
                                                                                                      ghost_fluid_extender_( material_manager_, halo_manager ),
                                                                                                      interface_extender_( halo_manager ) {
      /** Empty Constructor, besides initializer list. */
   }

public:
   MultiPhaseManager()                           = delete;
   ~MultiPhaseManager()                          = default;
   MultiPhaseManager( MultiPhaseManager const& ) = delete;
   MultiPhaseManager& operator=( MultiPhaseManager const& ) = delete;
   MultiPhaseManager( MultiPhaseManager&& )                 = delete;
   MultiPhaseManager& operator=( MultiPhaseManager&& ) = delete;

   /**
    * @brief Performs cut-cell mixing as provided by the CUT_CELL_MIXER object.
    * @param nodes The nodes which has to be mixed.
    */
   void Mix( std::vector<std::reference_wrapper<Node>> const& nodes ) const {
      static_cast<DerivedMultiPhaseManager const&>( *this ).MixImplementation( nodes );
   }

   /**
    * @brief Ensures a well-resolved level-set field satisfying the distance property.
    * @param nodes The node whose levelset field is considered.
    * @param is_last_stage Indicates whether scale separation has to be done or not (default = false).
    */
   void EnforceWellResolvedDistanceFunction( std::vector<std::reference_wrapper<Node>> const& nodes, bool const is_last_stage = false ) const {
      static_cast<DerivedMultiPhaseManager const&>( *this ).EnforceWellResolvedDistanceFunctionImplementation( nodes, is_last_stage );
   }

   /**
    * @brief Extend material states to ghost material.
    * @param nodes The nodes for which extension has to be done.
    */
   void Extend( std::vector<std::reference_wrapper<Node>> const& nodes ) const {
      static_cast<DerivedMultiPhaseManager const&>( *this ).ExtendPrimeStatesImplementation( nodes );
   }

   /**
    * @brief Extend interface quantities into narrow band.
    * @param nodes The nodes for which extension has to be done.
    */
   void ExtendInterfaceStates( std::vector<std::reference_wrapper<Node>> const& nodes ) const {
      static_cast<DerivedMultiPhaseManager const&>( *this ).ExtendInterfaceStatesImplementation( nodes );
   }

   /**
    * @brief After integration of the level-set field related quantities have to be adjusted to the propagated level-set field.
    *        Quantities which are adjusted are the interface tags and the volume fraction.
    * @param nodes The nodes on the finest level for which interface-tags and volume fractions have to be adjusted.
    * @param is_last_stage  Bool indicating whether this is the last stage of the RK cycle.
    */
   void UpdateIntegratedBuffer( std::vector<std::reference_wrapper<Node>> const& nodes, bool const is_last_stage ) const {
      static_cast<DerivedMultiPhaseManager const&>( *this ).UpdateIntegratedBufferImplementation( nodes, is_last_stage );
   }

   /**
    * @brief After integration of the level-set field related quantities have to be adjusted to the propagated level-set field.
    *        Quantities which are adjusted are the interface tags and the volume fraction.
    * @param nodes The nodes on the finest level for which interface-tags and volume fractions have to be adjusted.
    */
   void PropagateLevelset( std::vector<std::reference_wrapper<Node>> const& nodes ) const {
      static_cast<DerivedMultiPhaseManager const&>( *this ).PropagateLevelsetImplementation( nodes );
   }

   /**
    * @brief Computes interface quantities (as e.g. interface pressure and velocity) and extends them if necessary.
    * @param nodes A reference to the list containing the local multi nodes.
    * @param reset_interface_states Indicates whether interface quantities are set to zero before calculating an extending them.
    * @note Default value is false.
    */
   void ObtainInterfaceStates( std::vector<std::reference_wrapper<Node>> const& nodes, bool const reset_interface_states = false ) const {
      static_cast<DerivedMultiPhaseManager const&>( *this ).ObtainInterfaceStatesImplementation( nodes, reset_interface_states );
   }

   /**
    * @brief Initializes the volume fraction buffer. Volume fractions are calculated based on the reinitialized level-set buffer.
    * @note ATTENTION: This method should only be called during the initialization of the simulation.
    * @param nodes The nodes for which the volume fractions should be calculated.
    */
   void InitializeVolumeFractionBuffer( std::vector<std::reference_wrapper<Node>> const& nodes ) const {
      static_cast<DerivedMultiPhaseManager const&>( *this ).InitializeVolumeFractionBufferImplementation( nodes );
   }

   /**
    * @brief Transform given volume averaged conservatives to conservatives. This is done by a multiplication with the volume fraction.
    * @param node The node for which conservatives are calculated.
    */
   void TransformToConservatives( Node& node ) const {
      buffer_handler_.TransformToConservatives( node );
   }
};

#endif//MULTI_PHASE_MANAGER_H
