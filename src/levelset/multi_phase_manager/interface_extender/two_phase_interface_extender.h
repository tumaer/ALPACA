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
#ifndef TWO_PHASE_INTERFACE_EXTENDER_H
#define TWO_PHASE_INTERFACE_EXTENDER_H

#include "interface_extender.h"

template<InterfaceFieldType field_type>
class TwoPhaseInterfaceExtender : public InterfaceExtender<TwoPhaseInterfaceExtender<field_type>, field_type> {
   // Definition for shortening
   using InterfaceExtenderSpecification = InterfaceExtender<TwoPhaseInterfaceExtender<field_type>, field_type>;

   friend InterfaceExtenderSpecification;

private:
   // use some variables from base class
   using InterfaceExtenderSpecification::epsilon_;
   using InterfaceExtenderSpecification::field_type_;
   using InterfaceExtenderSpecification::number_of_convergence_tracking_quantities_;
   static constexpr unsigned int stencil_width_ = 2;
   static constexpr unsigned int i_offset_      = CC::HS() - stencil_width_;
   static constexpr unsigned int j_offset_      = CC::DIM() != Dimension::One ? CC::HS() - stencil_width_ : 0;
   static constexpr unsigned int k_offset_      = CC::DIM() == Dimension::Three ? CC::HS() - stencil_width_ : 0;
   static constexpr unsigned int repetition_    = CC::HS() / stencil_width_;

   /**
    * @brief Actual implementation of the iterative extension equation to extend interface fields into narrow band.
    * @param node The node which contains the phase which is extended.
    * @param convergence_tracking_quantities A vector holding information about the convergence status of the iterative extension method.
    */
   void IterativeExtension( Node& node, std::vector<double>& convergence_tracking_quantities ) const {

      double const( &levelset_reinitialized )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer( InterfaceDescription::Levelset );
      std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()]    = node.GetInterfaceTags();

      std::array<double, DTI( CC::DIM() )> rhs_contributions;

      double interface_field_change[CC::TCX()][CC::TCY()][CC::TCZ()];

      for( unsigned int field_index = 0; field_index < IF::NOFTE( field_type_ ); ++field_index ) {
         double const one_normalization_constant                     = 1.0 / ( convergence_tracking_quantities[field_index] + epsilon_ );
         double( &interface_field )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetFieldBuffer( field_type_, IF::FITE( field_type_, field_index ) );

         for( unsigned int iteration = 0; iteration < repetition_; ++iteration ) {
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                     interface_field_change[i][j][k] = 0.0;
                  }//k
               }   //j
            }      //i

            //Loop through all narrow-band cells which are no cut-cells to compute increment for iterative reinitialization
            for( unsigned int i = CC::FICX() - i_offset_; i <= CC::LICX() + i_offset_; ++i ) {
               for( unsigned int j = CC::FICY() - j_offset_; j <= CC::LICY() + j_offset_; ++j ) {
                  for( unsigned int k = CC::FICZ() - k_offset_; k <= CC::LICZ() + k_offset_; ++k ) {
                     /**
                       * We also extend in the reinitialization band in order to have better convergence behaviour (Reduce influence of implicitly imposed boundary
                       * conditions at the end of the narrow band). The convergence criteria is only checked for the extension band.
                       */
                     if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::ReinitializationBand ) && std::abs( interface_tags[i][j][k] ) > ITTI( IT::NewCutCell ) ) {

                        // compute normal
                        std::array<double, 3> const normal = GetNormal( levelset_reinitialized, i, j, k );

                        double const ls_sign = Signum( levelset_reinitialized[i][j][k] );

                        // compute extension terms and store in vector
                        // i-direction
                        // direction of the gradient from interface towards the material in both directions
                        int const n_i        = ( normal[0] * ls_sign < 0.0 ) - ( normal[0] * ls_sign > 0.0 );
                        rhs_contributions[0] = ( interface_field[i + n_i][j][k] - interface_field[i][j][k] );
                        rhs_contributions[0] *= std::abs( normal[0] );

                        // j-direction
                        if( CC::DIM() != Dimension::One ) {
                           int const n_j        = ( normal[1] * ls_sign < 0.0 ) - ( normal[1] * ls_sign > 0.0 );
                           rhs_contributions[1] = ( interface_field[i][j + n_j][k] - interface_field[i][j][k] );
                           rhs_contributions[1] *= std::abs( normal[1] );
                        }

                        // k-direction
                        if( CC::DIM() == Dimension::Three ) {
                           int const n_k        = ( normal[2] * ls_sign < 0.0 ) - ( normal[2] * ls_sign > 0.0 );
                           rhs_contributions[2] = ( interface_field[i][j][k + n_k] - interface_field[i][j][k] );
                           rhs_contributions[2] *= std::abs( normal[2] );
                        }

                        interface_field_change[i][j][k] = ConsistencyManagedSum( rhs_contributions ) * InterfaceStateExtensionConstants::Dtau;
                        if( InterfaceStateExtensionConstants::TrackConvergence && std::abs( interface_tags[i][j][k] ) <= ITTI( IT::ExtensionBand ) ) {
                           convergence_tracking_quantities[IF::NOFTE( field_type_ )] = std::max( convergence_tracking_quantities[IF::NOFTE( field_type_ )], std::abs( interface_field_change[i][j][k] * one_normalization_constant ) );
                        }
                     }//if cut cell
                  }   //k
               }      //j
            }         //i

            for( unsigned int i = CC::FICX() - i_offset_; i <= CC::LICX() + i_offset_; ++i ) {
               for( unsigned int j = CC::FICY() - j_offset_; j <= CC::LICY() + j_offset_; ++j ) {
                  for( unsigned int k = CC::FICZ() - k_offset_; k <= CC::LICZ() + k_offset_; ++k ) {
                     interface_field[i][j][k] += interface_field_change[i][j][k];
                  }//k
               }   //j
            }      //i
         }

      }//field of interface field_type
   }

public:
   TwoPhaseInterfaceExtender() = delete;
   explicit TwoPhaseInterfaceExtender( HaloManager& halo_manager ) : InterfaceExtenderSpecification( halo_manager ) {
      /** Empty besides initilializer list for call of base class constructor */
   }
   ~TwoPhaseInterfaceExtender()                                  = default;
   TwoPhaseInterfaceExtender( TwoPhaseInterfaceExtender const& ) = delete;
   TwoPhaseInterfaceExtender& operator=( TwoPhaseInterfaceExtender const& ) = delete;
   TwoPhaseInterfaceExtender( TwoPhaseInterfaceExtender&& )                 = delete;
   TwoPhaseInterfaceExtender& operator=( TwoPhaseInterfaceExtender&& ) = delete;
};

#endif//TWO_PHASE_INTERFACE_EXTENDER_H
