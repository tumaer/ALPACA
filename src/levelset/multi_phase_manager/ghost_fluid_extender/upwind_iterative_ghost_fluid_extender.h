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
#ifndef UPWIND_GHOST_FLUID_EXTENDER_H
#define UPWIND_GHOST_FLUID_EXTENDER_H

#include "ghost_fluid_extender.h"

/**
 * @brief The UpwindGhostFluidExtender class extends material-states from real-material to ghost-material by iteratively solving the extension equation.
 *
 * @tparam field_type The Fluid field type for which the extender is used (PrimeStates, Parameters or Conservatives)
 */
template<MaterialFieldType field_type>
class UpwindGhostFluidExtender : public GhostFluidExtender<UpwindGhostFluidExtender<field_type>, field_type> {
   // Definition for shortening
   using GhostFluidExtenderSpecification = GhostFluidExtender<UpwindGhostFluidExtender<field_type>, field_type>;
   // friend declaration
   friend GhostFluidExtenderSpecification;

private:
   // Taking some variables from the base class
   using GhostFluidExtenderSpecification::epsilon_;
   using GhostFluidExtenderSpecification::field_type_;
   using GhostFluidExtenderSpecification::number_of_convergence_tracking_quantities_;
   static constexpr unsigned int stencil_width_ = 2;
   static constexpr unsigned int i_offset_      = CC::HS() - stencil_width_;
   static constexpr unsigned int j_offset_      = CC::DIM() != Dimension::One ? CC::HS() - stencil_width_ : 0;
   static constexpr unsigned int k_offset_      = CC::DIM() == Dimension::Three ? CC::HS() - stencil_width_ : 0;
   static constexpr unsigned int repetition_    = CC::HS() / stencil_width_;

protected:
   /**
    * @brief Extend to cut-cell neighbors and extension band iteratively.
    * @param node The node which is extended.
    * @param convergence_tracking_quantities An array holding information about the convergence status of the iterative extension method.
    */
   void IterativeExtension( Node& node, double ( &convergence_tracking_quantities )[2][number_of_convergence_tracking_quantities_] ) const {

      std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
      double const( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()]            = node.GetInterfaceBlock().GetReinitializedBuffer( InterfaceDescription::Levelset );
      double const( &volume_fraction )[CC::TCX()][CC::TCY()][CC::TCZ()]     = node.GetInterfaceBlock().GetBaseBuffer( InterfaceDescription::VolumeFraction );

      // Loop through all materials of the node
      for( auto& phase : node.GetPhases() ) {

         std::int8_t const material_sign        = MaterialSignCapsule::SignOfMaterial( phase.first );
         unsigned int const material_index      = phase.first == MaterialSignCapsule::PositiveMaterial() ? 0 : 1;
         double const reference_volume_fraction = ( material_sign > 0 ) ? 0.0 : 1.0;
         double const material_sign_double      = double( material_sign );

         double extension_rhs[MF::ANOF( field_type_ )][CC::TCX()][CC::TCY()][CC::TCZ()];
         double one_normalization_constant[MF::ANOF( field_type_ )];
         for( unsigned int field_index = 0; field_index < MF::ANOF( field_type_ ); ++field_index ) {
            one_normalization_constant[field_index] = 1.0 / ( std::max( epsilon_, convergence_tracking_quantities[material_index][field_index] ) );
         }

         for( unsigned int iteration = 0; iteration < repetition_; ++iteration ) {
            /**
             * Setting the extension_rhs buffer to zero is crucial!
             */
            for( unsigned int field_index = 0; field_index < MF::ANOF( field_type_ ); field_index++ ) {
               for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                  for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                     for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        extension_rhs[field_index][i][j][k] = 0.0;
                     }// k
                  }   // j
               }      // i
            }         // fields of field_type

            std::array<double, DTI( CC::DIM() )> rhs_contributions;

            unsigned int derivative_indices[DTI( CC::DIM() )][2];
            for( unsigned int d = 0; d < DTI( CC::DIM() ); ++d ) {
               for( unsigned int l = 0; l < 2; ++l ) {
                  derivative_indices[d][l] = 0;
               }
            }

            // Loop through internal block - finally, fill extension band and cut-cell neighbors
            for( unsigned int i = CC::FICX() - i_offset_; i <= CC::LICX() + i_offset_; ++i ) {
               for( unsigned int j = CC::FICY() - j_offset_; j <= CC::LICY() + j_offset_; ++j ) {
                  for( unsigned int k = CC::FICZ() - k_offset_; k <= CC::LICZ() + k_offset_; ++k ) {
                     double const cell_volume_fraction = reference_volume_fraction + material_sign_double * volume_fraction[i][j][k];
                     /**
                      * We also extend in the reinitialization band in order to have better convergence behaviour (Reduce influence of implicitly imposed boundary
                      * conditions at the end of the narrow band). The convergence criteria is only checked for the extension band.
                      */
                     if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::ReinitializationBand ) && ( cell_volume_fraction <= CC::ETH() ) ) {

                        std::array<double, 3> const normal = GetNormal( levelset, i, j, k, material_sign );

                        if( normal[0] > 0.0 ) {
                           derivative_indices[0][0] = i + 1;
                           derivative_indices[0][1] = i;
                        } else {
                           derivative_indices[0][0] = i;
                           derivative_indices[0][1] = i - 1;
                        }

                        if constexpr( CC::DIM() != Dimension::One ) {
                           if( normal[1] > 0.0 ) {
                              derivative_indices[1][0] = j + 1;
                              derivative_indices[1][1] = j;
                           } else {
                              derivative_indices[1][0] = j;
                              derivative_indices[1][1] = j - 1;
                           }
                        }

                        if constexpr( CC::DIM() == Dimension::Three ) {
                           if( normal[2] > 0.0 ) {
                              derivative_indices[2][0] = k + 1;
                              derivative_indices[2][1] = k;
                           } else {
                              derivative_indices[2][0] = k;
                              derivative_indices[2][1] = k - 1;
                           }
                        }

                        //calculate gradients
                        for( unsigned int field_index = 0; field_index < MF::ANOF( field_type_ ); field_index++ ) {
                           double const( &cell )[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.second.GetFieldBuffer( field_type_, field_index );

                           rhs_contributions[0] = cell[derivative_indices[0][0]][j][k] - cell[derivative_indices[0][1]][j][k];
                           rhs_contributions[0] *= levelset[derivative_indices[0][0]][j][k] - levelset[derivative_indices[0][1]][j][k];

                           if constexpr( CC::DIM() != Dimension::One ) {
                              rhs_contributions[1] = cell[i][derivative_indices[1][0]][k] - cell[i][derivative_indices[1][1]][k];
                              rhs_contributions[1] *= levelset[i][derivative_indices[1][0]][k] - levelset[i][derivative_indices[1][1]][k];
                           }

                           if constexpr( CC::DIM() == Dimension::Three ) {
                              rhs_contributions[2] = cell[i][j][derivative_indices[2][0]] - cell[i][j][derivative_indices[2][1]];
                              rhs_contributions[2] *= levelset[i][j][derivative_indices[2][0]] - levelset[i][j][derivative_indices[2][1]];
                           }

                           extension_rhs[field_index][i][j][k] = ConsistencyManagedSum( rhs_contributions ) * ExtensionConstants::Dtau * material_sign_double;

                           if( ExtensionConstants::TrackConvergence && std::abs( interface_tags[i][j][k] ) <= ITTI( IT::ExtensionBand ) ) {
                              convergence_tracking_quantities[material_index][MF::ANOF( field_type_ )] = std::max( convergence_tracking_quantities[material_index][MF::ANOF( field_type_ )], std::abs( extension_rhs[field_index][i][j][k] * one_normalization_constant[field_index] ) );
                           }

                        }// fields of field_type
                     }   // cells to extend
                  }      // k
               }         // j
            }            // i

            for( unsigned int field_index = 0; field_index < MF::ANOF( field_type_ ); field_index++ ) {
               double( &extension_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.second.GetFieldBuffer( field_type_, field_index );
               for( unsigned int i = CC::FICX() - i_offset_; i <= CC::LICX() + i_offset_; ++i ) {
                  for( unsigned int j = CC::FICY() - j_offset_; j <= CC::LICY() + j_offset_; ++j ) {
                     for( unsigned int k = CC::FICZ() - k_offset_; k <= CC::LICZ() + k_offset_; ++k ) {
                        extension_buffer[i][j][k] += extension_rhs[field_index][i][j][k];
                     }// k
                  }   // j
               }      // i
            }         // fields of field_type
         }
      }
   }

public:
   UpwindGhostFluidExtender() = delete;
   explicit UpwindGhostFluidExtender( MaterialManager const& material_manager, HaloManager& halo_manager ) : GhostFluidExtenderSpecification( material_manager, halo_manager ) {
      // Empty Constructor, besides call of base class constructor.
   }
   ~UpwindGhostFluidExtender()                                 = default;
   UpwindGhostFluidExtender( UpwindGhostFluidExtender const& ) = delete;
   UpwindGhostFluidExtender& operator=( UpwindGhostFluidExtender const& ) = delete;
   UpwindGhostFluidExtender( UpwindGhostFluidExtender&& )                 = delete;
   UpwindGhostFluidExtender operator=( UpwindGhostFluidExtender&& ) = delete;
};

#endif//UPWIND_GHOST_FLUID_EXTENDER_H
