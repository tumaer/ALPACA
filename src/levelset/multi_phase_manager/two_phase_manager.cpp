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
#include "two_phase_manager.h"

#include "interface_tags/interface_tag_functions.h"
#include "material_sign_capsule.h"
#include "utilities/buffer_operations_interface.h"

/**
 * @brief Default constructor for the TwoPhaseManager. Calls the default constructor of the base class.
 * @param setup Instance to get access to user defined properties, relevant for the simulation.
 * @param material_manager Instance of a material manager, which already has been initialized according to the user input.
 * @param halo_manager Instance to a HaloManager which provides MPI-related methods.
 */
TwoPhaseManager::TwoPhaseManager( MaterialManager const& material_manager, HaloManager& halo_manager ) : MultiPhaseManager( material_manager, halo_manager ) {
#ifndef PERFORMANCE
   // NH I think != 2 would block degeneration to single phase, but I'm not sure.
   if( material_manager.GetNumberOfMaterials() > 2 ) { throw std::logic_error( "Do not use TwoPhaseManager for more than two materials" ); }
#endif
}

/**
 * @brief Allow cut-cell mixing for a single level-set field satisfying the signed-distance property. See also base class.
 * @param nodes See base class.
 */
void TwoPhaseManager::MixImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const {

   // TODO-19 JW: If integration is done on total cells, this halo update can possibly be left out
   halo_manager_.MaterialHaloUpdateOnLmax( MaterialFieldType::Conservatives );

   for( Node& node : nodes ) {
      cut_cell_mixer_.Mix( node );
      buffer_handler_.TransformToVolumeAveragedConservatives( node );
   }

   halo_manager_.MaterialHaloUpdateOnLmax( MaterialFieldType::Conservatives );

   /** Extension is done on the prime state buffer. Thus, cells in which we do not extend have to be filled with prime states obtained from
    * the integrated conservatives from the right-hand side buffer. This leads to the following occupation of the prime-state buffers:
    * Real-material cells and cut-cells in which we do not extend: Prime states obtained from integrated conservatives.
    * Extension-band cells and cut-cells in which we extend: Prime states of the last RK stage.
    */
   for( Node& node : nodes ) {
      buffer_handler_.CalculatePrimesFromIntegratedConservatives( node );
   }
}

/**
 * @brief Ensures a well-resolved single-level set field satisfying the signed-distance property. Therefore, scale-separation and reinitialization are performed.
 * @param nodes See base class.
 * @param is_last_stage See base class.
 */
void TwoPhaseManager::EnforceWellResolvedDistanceFunctionImplementation( std::vector<std::reference_wrapper<Node>> const& nodes, bool const is_last_stage ) const {
   if( CC::ScaleSeparationActive() && is_last_stage ) {
      scale_separator_.SeparateScales( nodes );
      /**
       * Since we also want to reinitialize scale-separated cells we cannot update the interface tags at this place.
       * We have to do it after the reinitialization.
       */
   }

   bool reinitialize = true;
   if( ReinitializationConstants::ReinitializeOnlyInLastRkStage && !is_last_stage ) {
      reinitialize = false;
   }

   if( reinitialize ) {
      levelset_reinitializer_.Reinitialize( nodes, is_last_stage );
   }

   if( is_last_stage ) {
      for( Node& node : nodes ) {
         buffer_handler_.AdaptConservativesToWellResolvedDistanceFunction( node );
      }
      halo_manager_.MaterialHaloUpdateOnLmax( MaterialFieldType::Conservatives );
      UpdateInterfaceTagsOnFinestLevel( nodes );

      for( Node& node : nodes ) {
         SetVolumeFractionBuffer( node );
      }
      halo_manager_.InterfaceHaloUpdateOnLmax( InterfaceBlockBufferType::VolumeFractionBase );
   }
}

/**
 * @brief Allows to extend material states to ghost cells.
 * @param node See base class.
 */
void TwoPhaseManager::ExtendPrimeStatesImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const {

   /** After the extension we have the following occupation of the prime-state buffer:
    * Real-material cells and cut-cells in which we do not extend: Prime states obtained from integrated conservatives.
    * Extension-band cells and cut-cells in which we extend: Extended values.
    * Other cells (Non-real material or narrow band cells): 0.0. */
   ghost_fluid_extender_.Extend( nodes );

   /** After the extension the right-hand side buffer has to be populated with the extended values. This leads to the following occupation of the
    * right-hand side buffer:
    * Real-fluid cells and cut-cells in which we do not extend: Integrated conservatives.
    * Extension-band cells and cut-cells in which we extend: Conservatives obtained from the extended prime states.
    * Other cells (Non-real material or narrow band cells): 0.0. */
   for( Node& node : nodes ) {
      buffer_handler_.CalculateConservativesFromExtendedPrimes( node );
   }//nodes

   /** To have integrated and extended conservatives also in the halo cells, we have to do a Halo Update on the right-hand side buffer. */
   halo_manager_.MaterialHaloUpdateOnLmax( MaterialFieldType::Conservatives );
}

/**
 * @brief Extends interface parameter into narrow band.
 * @param nodes See base class.
 */
void TwoPhaseManager::ExtendInterfaceStatesImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const {

   // Calls the extension method
   // halo update of the interface quantities to extend (extender acts on interface states)
   for( InterfaceState const is : IF::ISTE() ) {
      halo_manager_.InterfaceHaloUpdateOnLmax( MapInterfaceStateToInterfaceBlockBufferType( is ) );
   }
   // Extend all interface states
   interface_extender_.Extend( nodes );
}

/**
 * @brief Calculates the volume fractions for the positive material on the inner cells of a given node and saves them in the volume-fraction buffer.
 *        The volume fractions are calculated using the reinitialized levelset buffer.
 * @param node The node for which the volume fraction buffer is set.
 */
void TwoPhaseManager::SetVolumeFractionBuffer( Node& node ) const {

   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();

   InterfaceBlock& interface_block                             = node.GetInterfaceBlock();
   double const( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()]  = interface_block.GetReinitializedBuffer( InterfaceDescription::Levelset );
   double( &volume_fraction )[CC::TCX()][CC::TCY()][CC::TCZ()] = interface_block.GetBaseBuffer( InterfaceDescription::VolumeFraction );

   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell ) ) {// volume fraction calculation only makes sense for old or new cut cells, for non-cut cells it is trivial.
               volume_fraction[i][j][k] = geometry_calculator_.ComputeVolumeFraction( levelset, i, j, k );
            } else if( interface_tags[i][j][k] > ITTI( IT::NewCutCell ) ) {
               volume_fraction[i][j][k] = 1.0;
            } else {
               volume_fraction[i][j][k] = 0.0;
            }
         }//k
      }   //j
   }      //i
}

/**
 * See base class.
 * @param nodes See base class.
 */
void TwoPhaseManager::PropagateLevelsetImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const {

   /***
    * At this point in the algorithm time integration of the level-set field is already done. To fully propagate the level-set field,
    * interface tags and volume fractions have to be updated.
    * The interface tag update and the calculation of volume fractions are based on the level-set values in the reinitialized level-set buffer. Thus, after the advection of the level-set field,
    * the advected level-set values, which are currently in the right-hand side level-set buffer, have to be copied to the reinitialized level-set buffer.
    */
   BO::Interface::CopyInterfaceDescriptionBufferForNodeList<InterfaceDescriptionBufferType::RightHandSide, InterfaceDescriptionBufferType::Reinitialized>( nodes );
   UpdateInterfaceTagsOnFinestLevel( nodes );

   // Set the volume fraction buffer according to the propagated level-set field.
   for( Node& node : nodes ) {
      SetVolumeFractionBuffer( node );
   }
   // A halo update for the volume fractions is necessary to have also correct volume fractions in the halo cells.
   halo_manager_.InterfaceHaloUpdateOnLmax( InterfaceBlockBufferType::VolumeFractionBase );
}

/**
 * @brief See base class.
 * @param nodes See base class.
 */
void TwoPhaseManager::InitializeVolumeFractionBufferImplementation( std::vector<std::reference_wrapper<Node>> const& nodes ) const {
   for( Node& node : nodes ) {
      SetVolumeFractionBuffer( node );
   }
   halo_manager_.InterfaceHaloUpdateOnLmax( InterfaceBlockBufferType::VolumeFractionBase );
}

/**
 * @brief Sets the interface tags on the finest level.
 * @param nodes The nodes on the finest level, which have a level-set block.
 */
void TwoPhaseManager::UpdateInterfaceTagsOnFinestLevel( std::vector<std::reference_wrapper<Node>> const& nodes ) const {

   for( Node& node : nodes ) {
      InterfaceTagFunctions::SetInternalCutCellTagsFromLevelset( node.GetInterfaceBlock().GetReinitializedBuffer( InterfaceDescription::Levelset ), node.GetInterfaceTags() );
   }
   halo_manager_.InterfaceTagHaloUpdateOnLmax();

   for( Node& node : nodes ) {
      InterfaceTagFunctions::SetTotalInterfaceTagsFromCutCells( node.GetInterfaceTags() );
   }
   halo_manager_.InterfaceTagHaloUpdateOnLmax();
}

/**
 * @brief See base class.
 * @param nodes See base class.
 */
void TwoPhaseManager::ObtainInterfaceStatesImplementation( std::vector<std::reference_wrapper<Node>> const& nodes, bool const reset_interface_states ) const {
   if( reset_interface_states ) {
      for( Node& node : nodes ) {
         BO::SetFieldBuffer( node.GetInterfaceBlock().GetInterfaceStateBuffer(), 0.0 );
      }
   }
   for( auto const& node : nodes ) {
      // solve interface Riemann problem to obtain interface velocity and interface exchange terms
      interface_state_calculator_.ObtainInterfaceStates( node );
   }

   for( InterfaceState const is : IF::ASOS() ) {
      halo_manager_.InterfaceHaloUpdateOnLmax( MapInterfaceStateToInterfaceBlockBufferType( is ) );
   }

   if constexpr( InterfaceStateTreatmentConstants::ExtendInterfaceQuantities ) {
      //extends interface quantities
      ExtendInterfaceStates( nodes );
   }
}