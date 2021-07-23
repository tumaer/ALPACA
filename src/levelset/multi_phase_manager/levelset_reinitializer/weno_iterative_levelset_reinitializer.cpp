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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#include "weno_iterative_levelset_reinitializer.h"
#include "user_specifications/two_phase_constants.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief Computes the advection velocity used for the iterative reinitialization. The advection velocity is a smoothened signum of the levelset.
 *        Smoothening is necessary for small absolute level-set values to not reinitialize there.
 * @param levelset The levelset-value.
 * @return The advection velocity.
 */
inline double ComputeAdvectionVelocity( double const levelset ) {
   constexpr double epsilon = 1.0;
   return levelset / std::sqrt( levelset * levelset + ( epsilon * epsilon ) );
}

/**
 * @brief The default constructor for a WenoIterativeLevelsetReinitializer. Calls the default constructor of the base class.
 * @param halo_manager See base class.
 */
WenoIterativeLevelsetReinitializer::WenoIterativeLevelsetReinitializer( HaloManager& halo_manager ) : IterativeLevelsetReinitializerBase( halo_manager ) {
   // Empty Constructor, besides call of base class constructor.
}

/**
 * @brief Reinitializes the levelset field of a node using a HJ WENO scheme.
 * @param node The node with levelset block which has to be reinitialized.
 * @param levelset_type  Level-set field type that is reinitialized.
 * @param is_last_stage Return whether it's the last RK stage or not.
 * @return The residuum for the current node.
 */
double WenoIterativeLevelsetReinitializer::ReinitializeSingleNodeImplementation( Node& node, InterfaceDescriptionBufferType const levelset_type, bool const is_last_stage ) const {

   using ReconstructionStencil = ReconstructionStencilSetup::Concretize<levelset_reconstruction_stencil>::type;

   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags( levelset_type );
   InterfaceBlock& interface_block                                       = node.GetInterfaceBlock();
   double( &levelset_orig )[CC::TCX()][CC::TCY()][CC::TCZ()]             = interface_block.GetInterfaceDescriptionBuffer( levelset_type )[InterfaceDescription::Levelset];
   double const( &levelset_0_orig )[CC::TCX()][CC::TCY()][CC::TCZ()]     = interface_block.GetRightHandSideBuffer( InterfaceDescription::Levelset );

   double reinitialization_rhs[CC::TCX()][CC::TCY()][CC::TCZ()];

   double residuum = 0.0;

   double derivatives[DTI( CC::DIM() )][2];
   for( unsigned int d = 0; d < DTI( CC::DIM() ); ++d ) {
      for( unsigned int e = 0; e < 2; ++e ) {
         derivatives[d][e] = 0.0;
      }
   }

   for( unsigned int i = 0; i < CC::TCX(); ++i ) {
      for( unsigned int j = 0; j < CC::TCY(); ++j ) {
         for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
            reinitialization_rhs[i][j][k] = 0.0;
         }// k
      }   // j
   }      // i

   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            if( std::abs( interface_tags[i][j][k] ) > ITTI( IT::NewCutCell ) || ( ReinitializationConstants::ReinitializeCutCells && is_last_stage ) ) {

               // We normalize the level-set field on the cell size, thus 1.0 is necessary as cell size for stencil evaluation.
               derivatives[0][0] = SU::Derivative<ReconstructionStencil, SP::UpwindLeft, Direction::X>( levelset_orig, i, j, k, 1.0 );
               derivatives[0][1] = SU::Derivative<ReconstructionStencil, SP::UpwindRight, Direction::X>( levelset_orig, i, j, k, 1.0 );

               if constexpr( CC::DIM() != Dimension::One ) {
                  derivatives[1][0] = SU::Derivative<ReconstructionStencil, SP::UpwindLeft, Direction::Y>( levelset_orig, i, j, k, 1.0 );
                  derivatives[1][1] = SU::Derivative<ReconstructionStencil, SP::UpwindRight, Direction::Y>( levelset_orig, i, j, k, 1.0 );
               }

               if constexpr( CC::DIM() == Dimension::Three ) {
                  derivatives[2][0] = SU::Derivative<ReconstructionStencil, SP::UpwindLeft, Direction::Z>( levelset_orig, i, j, k, 1.0 );
                  derivatives[2][1] = SU::Derivative<ReconstructionStencil, SP::UpwindRight, Direction::Z>( levelset_orig, i, j, k, 1.0 );
               }

               double const old_levelset_sign   = Signum( levelset_0_orig[i][j][k] );
               double const godunov_hamiltonian = GodunovHamiltonian( derivatives, old_levelset_sign );
               double const advection_velocity  = ComputeAdvectionVelocity( levelset_0_orig[i][j][k] );
               double const increment           = ReinitializationConstants::Dtau * advection_velocity * ( 1.0 - godunov_hamiltonian );
               if( ReinitializationConstants::TrackConvergence && std::abs( interface_tags[i][j][k] ) < ITTI( IT::ReinitializationBand ) ) {
                  residuum = std::max( residuum, std::abs( increment ) );
               }

               reinitialization_rhs[i][j][k] = increment;
            }
         }//k
      }   //j
   }      //i

   //add to level-set and also the residuum
   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            levelset_orig[i][j][k] += reinitialization_rhs[i][j][k];
         }//k
      }   //j
   }      //i

   return residuum;
}
