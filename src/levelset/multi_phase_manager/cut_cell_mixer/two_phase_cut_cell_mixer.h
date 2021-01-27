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
#ifndef TWO_PHASE_CUT_CELL_MIXER_H
#define TWO_PHASE_CUT_CELL_MIXER_H

#include "cut_cell_mixer.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "enums/interface_tag_definition.h"

/**
 * @brief Class that implements the main functionality for mixing between two fluids. Specializations
 *        necessary since different choices for the mixing contributions exist.
 * @tparam DerivedTwoPhaseCutCellMixer Typename as template parameter due to CRTP.
 */
template<typename DerivedTwoPhaseCutCellMixer>
class TwoPhaseCutCellMixer : public CutCellMixer<DerivedTwoPhaseCutCellMixer> {

   friend CutCellMixer<DerivedTwoPhaseCutCellMixer>;
   friend DerivedTwoPhaseCutCellMixer;

private:
   unsigned int const number_of_mixing_operations_;
   using CutCellMixer<DerivedTwoPhaseCutCellMixer>::material_manager_;

   /**
     * @brief FICMO Gives the index of the First Internal Cell Minus One in a block per dimension.
     */
   static constexpr unsigned int FICMOX = CC::HS() - 1;
   static constexpr unsigned int FICMOY = CC::DIM() != Dimension::One ? CC::HS() - 1 : 0;
   static constexpr unsigned int FICMOZ = CC::DIM() == Dimension::Three ? CC::HS() - 1 : 0;

   /**
    * @brief LICPO Gives the index of the Last Internal Cell Plus One in a block per dimension. I. e. the returned index must be included if the internal cells are of interest.
    */
   static constexpr unsigned int LICPOX = CC::HS() + CC::ICX();
   static constexpr unsigned int LICPOY = CC::DIM() != Dimension::One ? ( CC::HS() + CC::ICY() ) : 0;
   static constexpr unsigned int LICPOZ = CC::DIM() == Dimension::Three ? ( CC::HS() + CC::ICZ() ) : 0;

   /**
    * @brief The default constructor for the TwoPhaseCutCellMixer class.
    * @param halo_manager Instance to a HaloManager which provides MPI-related methods.
    */
   explicit TwoPhaseCutCellMixer( HaloManager& halo_manager, unsigned int const number_of_mixing_operations, MaterialManager const& material_manager ) : CutCellMixer<DerivedTwoPhaseCutCellMixer>( halo_manager, material_manager ),
                                                                                                                                                         number_of_mixing_operations_( number_of_mixing_operations ) {
      // Empty Constructor, besides call of base class constructor.
   }

   /**
    * @brief Performs mixing on given nodes. Mixing weights are obtained by a interface-normal based splitting.
    * @param node See base class.
    */
   void MixImplementation( Node& node ) const {

      /***
      * mixing_contributions contains the relevant information about the mixing operations of one block.
      * The outermost vector separates the mixing contributions cell-wise.
      * The following pair is to separate unsigned-int and double-valued information about a mixing operations.
      * Each element of the vector of a pair contains information about a single mixing operation.
      * First vector (unsigned int): i, j, k, i_target, j_target, k_target.
      * Second vector (double): mixing fraction, factor to calculate mixing fluxes.
      */
      std::vector<std::pair<std::vector<std::array<unsigned int, 6>>, std::vector<std::array<double, 2>>>> mixing_contributions;
      double mixing_fluxes[MF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()];

      for( auto const& phase : node.GetPhases() ) {
         MaterialName const material = phase.first;
         if constexpr( CC::SolidBoundaryActive() ) {
            if( material_manager_.IsSolidBoundary( material ) ) {
               continue;
            }
         }
         //reset mixing fluxes and mixing contributions for the respective phase
         for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                     mixing_fluxes[e][i][j][k] = 0.0;
                  }//k
               }   //j
            }      //i
         }         //equation
         mixing_contributions.clear();

         static_cast<DerivedTwoPhaseCutCellMixer const&>( *this ).CalculateMixingContributionsImplementation( node, material, mixing_contributions );
         for( unsigned int c = 0; c < number_of_mixing_operations_; ++c ) {
            CalculateMixingFluxes( node, material, mixing_contributions, c, mixing_fluxes );
         }
         AddMixingFluxesToConservatives( node, material, mixing_fluxes );
      }//phases
   }

   /**
    * @brief Adds the mixing fluxes saved in conservative_change to the right-hand side conservative buffer of a specific phase.
    * @param node The node which contains the phase that has to be mixed.
    * @param material The material which specifies the phase that has to be mixed.
    * @param mixing_fluxes The mixing fluxes which are added to the right-hand side buffer.
    */
   void AddMixingFluxesToConservatives( Node& node, MaterialName const material, double const ( &mixing_fluxes )[MF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()] ) const {

      Block& phase                                                          = node.GetPhaseByMaterial( material );
      std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

      for( Equation const eq : MF::ASOE() ) {
         double( &conservatives )[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.GetRightHandSideBuffer( eq );
         for( unsigned int i = FICMOX; i <= LICPOX; ++i ) {
            for( unsigned int j = FICMOY; j <= LICPOY; ++j ) {
               for( unsigned int k = FICMOZ; k <= LICPOZ; ++k ) {
                  if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::ExtensionBand ) ) {
                     conservatives[i][j][k] += mixing_fluxes[ETI( eq )][i][j][k];
                  }
               }//k
            }   //j
         }      //i
      }         //equation
   }

   /**
    * @brief Based on given mixing contributions, i.e. cell-pairs which are mixed, the mixing fluxes are calculated.
    * @param node The node for which the mixing fluxes are calculated.
    * @param material The material which specifies the phase for which the mixing contributions are calculated.
    * @param mixing_contributions The mixing contributions for which the corresponding fluxes are calculated. A mixing contribution tuple contains information about the mixing operation
    *        between two cells. It contains the indices of the target cell (i_target, j_target and k_target - saved as unsigned int), the mixing fraction beta (saved as double) and a factor necessary
    *        to calculate the mixing fluxes (saved a a double).
    * @param contribution_identifier An identifier which mixing contribution is considered.
    * @param mixing_fluxes Indirect return parameter for the mixing fluxes.
    */
   void CalculateMixingFluxes( Node const& node, MaterialName const material,
                               std::vector<std::pair<std::vector<std::array<unsigned int, 6>>, std::vector<std::array<double, 2>>>> const& mixing_contributions,
                               unsigned int const contribution_identifier,
                               double ( &mixing_fluxes )[MF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()] ) const {

      Block const& phase                                                = node.GetPhaseByMaterial( material );
      std::int8_t const material_sign                                   = MaterialSignCapsule::SignOfMaterial( material );
      double const( &volume_fraction )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer( InterfaceDescription::VolumeFraction );

      double const reference_volume_fraction = ( material_sign > 0 ) ? 0.0 : 1.0;
      double const material_sign_double      = double( material_sign );

      for( auto const& mix_element : mixing_contributions ) {
         std::vector<std::array<unsigned int, 6>> const& indices = mix_element.first;
         std::vector<std::array<double, 2>> const& flux_weights  = mix_element.second;
         if( contribution_identifier < indices.size() ) {
            unsigned int const i        = indices[contribution_identifier][0];
            unsigned int const j        = indices[contribution_identifier][1];
            unsigned int const k        = indices[contribution_identifier][2];
            unsigned int const i_target = indices[contribution_identifier][3];
            unsigned int const j_target = indices[contribution_identifier][4];
            unsigned int const k_target = indices[contribution_identifier][5];
            double const factor         = flux_weights[contribution_identifier][1];

            double const volume_fraction_target = reference_volume_fraction + material_sign_double * volume_fraction[i_target][j_target][k_target];
            double const volume_fraction_self   = reference_volume_fraction + material_sign_double * volume_fraction[i][j][k];

            for( Equation const eq : MF::ASOE() ) {
               double const conservative_target = phase.GetRightHandSideBuffer( eq )[i_target][j_target][k_target];
               double const conservative_self   = phase.GetRightHandSideBuffer( eq )[i][j][k];
               double const M                   = factor * ( conservative_target * volume_fraction_self - conservative_self * volume_fraction_target );

               mixing_fluxes[ETI( eq )][i][j][k] += M;
               mixing_fluxes[ETI( eq )][i_target][j_target][k_target] -= M;
            }
         }
      }
   }

public:
   TwoPhaseCutCellMixer()                              = delete;
   ~TwoPhaseCutCellMixer()                             = default;
   TwoPhaseCutCellMixer( TwoPhaseCutCellMixer const& ) = delete;
   TwoPhaseCutCellMixer& operator=( TwoPhaseCutCellMixer const& ) = delete;
   TwoPhaseCutCellMixer( TwoPhaseCutCellMixer&& )                 = delete;
   TwoPhaseCutCellMixer& operator=( TwoPhaseCutCellMixer&& ) = delete;
};

#endif//TWO_PHASE_CUT_CELL_MIXER_H
