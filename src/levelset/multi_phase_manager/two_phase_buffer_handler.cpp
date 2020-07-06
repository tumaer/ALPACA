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
#include "two_phase_buffer_handler.h"

#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "enums/interface_tag_definition.h"

/**
 * @brief The constructor of the TwoPhaseBufferHandler class.
 * @param material_manager An instance to the material manager.
 */
TwoPhaseBufferHandler::TwoPhaseBufferHandler( MaterialManager const& material_manager ) : BufferHandler( material_manager ),
   prime_state_handler_( material_manager )
{
   //Empty besides call of base class constructor
}

/**
 * @brief See base class.
 * @param node See base class.
 */
void TwoPhaseBufferHandler::TransformToConservativesImplementation( Node& node ) const {

   if( node.HasLevelset() ) {
      double const (&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetBaseBuffer( InterfaceDescription::VolumeFraction );
      for( auto& phase : node.GetPhases() ) {
         // for cut cells factorize volume fractions into conservatives -> volume averaged to real conservatives
         auto const material_sign = MaterialSignCapsule::SignOfMaterial( phase.first );
         double const reference_volume_fraction = ( material_sign > 0 ) ? 0.0 : 1.0;
         double const material_sign_double = double( material_sign );
         for( Equation const eq : MF::ASOE() ) {
            double (&average_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.second.GetAverageBuffer( eq );
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                     average_buffer[i][j][k] *= reference_volume_fraction + material_sign_double * volume_fraction[i][j][k];
                  } //k
               } //j
            } //i
         } //equation
      } //phases
   } //node contains levelset
}

/**
 * @brief See base class.
 * @param node See base class.
 */
void TwoPhaseBufferHandler::TransformToVolumeAveragedConservativesImplementation( Node& node ) const {
   InterfaceBlock const& interface_block = node.GetInterfaceBlock();
   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   double const (&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] = interface_block.GetBaseBuffer( InterfaceDescription::VolumeFraction );

   for( auto& phase : node.GetPhases() ) {
      std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial( phase.first );
      double const reference_volume_fraction = ( material_sign > 0 ) ? 0.0 : 1.0;
      double const material_sign_double = double( material_sign );
      for( Equation const eq : MF::ASOE() ) {
         double (&conservative)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.second.GetRightHandSideBuffer( eq );
         for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                  if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell) ) {
                     double const cell_volume_fraction = reference_volume_fraction + material_sign_double * volume_fraction[i][j][k];
                     if(cell_volume_fraction != 0.0) {
                        conservative[i][j][k] /= cell_volume_fraction;
                     } else {
                        conservative[i][j][k] = 0.0;
                     }
                  }
               } //k
            } //j
         } //i
      } //equation
   } //phases
}

/**
 * @brief See base class.
 * @param node See base class.
 */
void TwoPhaseBufferHandler::AdaptConservativesToWellResolvedDistanceFunctionImplementation( Node &node ) const {
   for( auto& phase : node.GetPhases() ) {
      std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
      double const (&volume_fraction_no_scale)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetBaseBuffer( InterfaceDescription::VolumeFraction );
      double const (&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer( InterfaceDescription::Levelset );

      MaterialName const material = phase.first;
      std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial( material );

      double const reference_volume_fraction = ( material_sign > 0 ) ? 0.0 : 1.0;
      double const material_sign_double = double( material_sign );

      Conservatives& conservatives_rhs = phase.second.GetRightHandSideBuffer();
      PrimeStates const& prime_states = phase.second.GetPrimeStateBuffer();

      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
               if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::ExtensionBand ) || std::abs( interface_tags[i][j][k] ) == ITTI( IT::ScaleSeparatedCell) ) {
                  double const alpha_no_scale = reference_volume_fraction + material_sign_double * volume_fraction_no_scale[i][j][k];
                  if( alpha_no_scale <= 0.0 && material_sign * levelset_reinitialized[i][j][k] > 0.0 ) { // TODO-19 JW: Check whether == 0.0 is suitable or <= epsilon should be used.
                     prime_state_handler_.ConvertPrimeStatesToConservatives( material, prime_states, conservatives_rhs, i, j, k );
                  }
               } //scale separated cells
            } //k
         } //j
      } //i
   } //phases
}

/**
 * @brief See base class.
 * @param node See base class.
 */
void TwoPhaseBufferHandler::CalculatePrimesFromIntegratedConservativesImplementation( Node &node ) const {
   InterfaceBlock const& interface_block = node.GetInterfaceBlock();
   double const (&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] = interface_block.GetBaseBuffer( InterfaceDescription::VolumeFraction );
   double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] = interface_block.GetReinitializedBuffer( InterfaceDescription::Levelset );
   for( auto& phase : node.GetPhases() ) {
      PrimeStates& prime_states = phase.second.GetPrimeStateBuffer();
      MaterialName const material = phase.first;
      std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial( material );
      double const reference_volume_fraction = ( material_sign > 0 ) ? 0.0 : 1.0;
      double const material_sign_double = double( material_sign );
      Conservatives const& conservatives_rhs = phase.second.GetRightHandSideBuffer();
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               double const cell_volume_fraction = reference_volume_fraction + material_sign_double * volume_fraction[i][j][k];
               if( material_sign * levelset[i][j][k] > 0.0 && cell_volume_fraction > CC::MITH() ) {
                  prime_state_handler_.ConvertConservativesToPrimeStates( material, conservatives_rhs, prime_states, i, j, k );
               } //cells in which is not extended
            } //k
         } //j
      } //i
   } //phases
}

/**
 * @brief See base class.
 * @param node See base class.
 */
void TwoPhaseBufferHandler::CalculateConservativesFromExtendedPrimesImplementation( Node &node ) const {
   InterfaceBlock const& interface_block = node.GetInterfaceBlock();
   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] = interface_block.GetReinitializedBuffer( InterfaceDescription::Levelset );
   for( auto& phase : node.GetPhases() ) {
      MaterialName const material = phase.first;
      Conservatives& conservatives_rhs = phase.second.GetRightHandSideBuffer();
      std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial( material );
      PrimeStates const& prime_states = phase.second.GetPrimeStateBuffer();
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
               if( std::abs( interface_tags[i][j][k] ) <= ITTI(IT::ExtensionBand) || material_sign * levelset[i][j][k] > 0.0 ) {
                  //narrow band cells, in which is extended
                  prime_state_handler_.ConvertPrimeStatesToConservatives( material, prime_states, conservatives_rhs, i, j, k );
               } else if( material_sign * levelset[i][j][k] < 0.0 ) {
                  //ghost fluid cells in which is not extended
                  for( Equation const e : MF::ASOE() ) {
                     conservatives_rhs[e][i][j][k] = 0.0;
                  }
               }
            } //k
         } //j
      } //i
   } //phases
}
