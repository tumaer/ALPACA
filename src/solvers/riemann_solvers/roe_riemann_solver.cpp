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
#include "solvers/riemann_solvers/roe_riemann_solver.h"
#include "stencils/stencil_utilities.h"

namespace {
   /**
 * @brief Computes the advection within the provided block.
 * @param block Block of the phase under consideration.
 * @param advection Reference to an array which will be filled with the advection (indirect return parameter).
 * @tparam DIR Indicates which spatial direction is to be computed.
 * @note Hotpath function.
 */
   template<Direction DIR>
   void ComputeAdvection( Block const& block, double ( &advection )[MF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()] ) {

      PrimeStates const& prime_states                                      = block.GetPrimeStateBuffer();
      double const( &energy )[CC::TCX()][CC::TCY()][CC::TCZ()]             = block.GetAverageBuffer( Equation::Energy );
      double const( &direction_velocity )[CC::TCX()][CC::TCY()][CC::TCZ()] = prime_states[MF::AV()[DTI( DIR )]];

      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               double const cell_rho      = prime_states[PrimeState::Density][i][j][k];
               double const cell_pressure = prime_states[PrimeState::Pressure][i][j][k];
               double const cell_energy   = energy[i][j][k];

               advection[ETI( Equation::Mass )][i][j][k]   = cell_rho * direction_velocity[i][j][k];
               advection[ETI( Equation::Energy )][i][j][k] = ( cell_energy + cell_pressure ) * direction_velocity[i][j][k];

               for( unsigned int d = 0; d < DTI( CC::DIM() ); ++d ) {
                  advection[ETI( MF::AME()[d] )][i][j][k] = cell_rho * ( direction_velocity[i][j][k] * prime_states[MF::AV()[d]][i][j][k] );
               }
               // Add pressure to convective flux-term of momentum
               advection[ETI( MF::AME()[DTI( DIR )] )][i][j][k] += cell_pressure;
            }//Z-Loop
         }   //Y-Loop
      }      //X-Loop
   }
}// namespace

/**
 * @brief Standard constructor using an already existing MaterialManager and EigenDecomposition object.
 * @param material_manager .
 * @param eigendecomposition_calculator .
 */
RoeRiemannSolver::RoeRiemannSolver( MaterialManager const& material_manager, EigenDecomposition const& eigendecomposition_calculator ) : RiemannSolver( material_manager, eigendecomposition_calculator ) {
   /* Empty besides initializer list*/
}

/**
 * @brief Solving the right hand side of the Euler Equations. Using Roe transformation and flux splitting
 *        with the set stencil. Also See base class.
 * @note Hotpath function.
 */
void RoeRiemannSolver::UpdateImplementation(
      std::pair<MaterialName const, Block> const& mat_block, double const cell_size,
      double ( &fluxes_x )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double ( &fluxes_y )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double ( &fluxes_z )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1] ) const {

   double advection_contribution[MF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()];

   double roe_eigenvectors_left[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()];
   double roe_eigenvectors_right[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()];
   double fluxfunction_wavespeeds[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()];

   // NH This gives performance. Not sure why, propably first touch 'problems'.
   for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               advection_contribution[e][i][j][k] = 0.0;
            }
         }
      }
   }

   // Setting the eigenvectors, eigenvalues to zero is necessary for two-phase simulations as not every entry is necessarily set during eigendecomposition
   for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
      for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
         for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
            for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
               for( unsigned int f = 0; f < MF::ANOE(); ++f ) {
                  roe_eigenvectors_left[i][j][k][e][f]  = 0.0;
                  roe_eigenvectors_right[i][j][k][e][f] = 0.0;
               }
               fluxfunction_wavespeeds[i][j][k][e] = 0.0;
            }
         }
      }
   }

   eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::X>( mat_block, roe_eigenvectors_left, roe_eigenvectors_right, fluxfunction_wavespeeds );
   ComputeAdvection<Direction::X>( mat_block.second, advection_contribution );
   ComputeFluxes<Direction::X>( mat_block.second, fluxes_x, advection_contribution, cell_size, roe_eigenvectors_left, roe_eigenvectors_right, fluxfunction_wavespeeds );

   if constexpr( CC::DIM() != Dimension::One ) {
      eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::Y>( mat_block, roe_eigenvectors_left, roe_eigenvectors_right, fluxfunction_wavespeeds );
      ComputeAdvection<Direction::Y>( mat_block.second, advection_contribution );
      ComputeFluxes<Direction::Y>( mat_block.second, fluxes_y, advection_contribution, cell_size, roe_eigenvectors_left, roe_eigenvectors_right, fluxfunction_wavespeeds );
   }

   if constexpr( CC::DIM() == Dimension::Three ) {
      eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::Z>( mat_block, roe_eigenvectors_left, roe_eigenvectors_right, fluxfunction_wavespeeds );
      ComputeAdvection<Direction::Z>( mat_block.second, advection_contribution );
      ComputeFluxes<Direction::Z>( mat_block.second, fluxes_z, advection_contribution, cell_size, roe_eigenvectors_left, roe_eigenvectors_right, fluxfunction_wavespeeds );
   }
}

/**
 * @brief Computes the cell-face fluxes with the selected stencil using a Roe-type finite-difference flux splitting.
 * @param block The block of the phase under consideration.
 * @param fluxes Reference to an array which is filled with the computed fluxes (indirect return parameter).
 * @param advection Reference to an array holding the advection in the current block.
 * @param cell_size .
 * @param roe_eigenvectors_left .
 * @param roe_eigenvectors_right .
 * @param roe_eigenvalues .
 * @tparam DIR Indicates which spatial direction is to be computed.
 * @note Hotpath function.
 */
template<Direction DIR>
void RoeRiemannSolver::ComputeFluxes( Block const& block,
                                      double ( &fluxes )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                                      double ( &advection )[MF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()], double const cell_size,
                                      double ( &roe_eigenvectors_left )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()],
                                      double ( &roe_eigenvectors_right )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()],
                                      double ( &fluxfunction_wavespeed )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()] ) const {

   using ReconstructionStencil = ReconstructionStencilSetup::Concretize<reconstruction_stencil>::type;

   // Index shift for eigenvalues and fluxes
   constexpr int offset_x = CC::FICX() - 1;
   constexpr int offset_y = CC::DIM() != Dimension::One ? CC::FICY() - 1 : -1;
   constexpr int offset_z = CC::DIM() == Dimension::Three ? CC::FICZ() - 1 : -1;

   constexpr unsigned int x_varying = DIR == Direction::X ? 1 : 0;
   constexpr unsigned int y_varying = DIR == Direction::Y ? 1 : 0;
   constexpr unsigned int z_varying = DIR == Direction::Z ? 1 : 0;

   unsigned int const reconstruction_stencil_downstream_size = ReconstructionStencil::DownstreamStencilSize();

   //NH Compiler likes loops counters to be fixed - so we help him.
   constexpr unsigned int x_start = DIR == Direction::X ? CC::FICX() - 1 : CC::FICX();
   constexpr unsigned int y_start = DIR == Direction::Y ? CC::FICY() - 1 : CC::FICY();
   constexpr unsigned int z_start = DIR == Direction::Z ? CC::FICZ() - 1 : CC::FICZ();
   constexpr unsigned int x_end   = CC::LICX();
   constexpr unsigned int y_end   = CC::LICY();
   constexpr unsigned int z_end   = CC::LICZ();

   std::array<double, ReconstructionStencil::StencilSize()> positive_characteristic_flux;
   std::array<double, ReconstructionStencil::StencilSize()> negative_characteristic_flux;

   std::array<double, MF::ANOE()> characteristic_flux;

   auto const& conservatives = block.GetAverageBuffer();

   for( unsigned int i = x_start; i <= x_end; ++i ) {
      for( unsigned int j = y_start; j <= y_end; ++j ) {
         for( unsigned int k = z_start; k <= z_end; ++k ) {
            // Shifted indices to match block index system and roe-eigenvalue index system
            int const i_index = i - offset_x;
            int const j_index = j - offset_y;
            int const k_index = k - offset_z;

            auto const& alpha_fluxvectorsplitting = fluxfunction_wavespeed[i_index][j_index][k_index];

            // reconstruct fluxes at face i+1/2 by characteristic decomposition and applying an WENO scheme to smoothen stencils
            for( unsigned int n = 0; n < MF::ANOE(); ++n ) {// n is index of characteristic field (eigenvalue, eigenvector)
               auto const& roe_eigenvector_left_temp = roe_eigenvectors_left[i_index][j_index][k_index][n];

               for( unsigned int m = 0; m < ReconstructionStencil::StencilSize(); ++m ) {
                  // This resetting is necessary!
                  positive_characteristic_flux[m] = 0.0;
                  negative_characteristic_flux[m] = 0.0;
                  for( unsigned int const l : conservative_equation_summation_sequence_[DTI( DIR )] ) {// l is index of conservative equation, iterated in symmetry-preserving sequence
                     // Compute characteristics for U and advection
                     double const u_characteristic_temp         = conservatives[l][i + x_varying * ( m - reconstruction_stencil_downstream_size )][j + y_varying * ( m - reconstruction_stencil_downstream_size )][k + z_varying * ( m - reconstruction_stencil_downstream_size )] * roe_eigenvector_left_temp[l];
                     double const advection_characteristic_temp = advection[l][i + x_varying * ( m - reconstruction_stencil_downstream_size )][j + y_varying * ( m - reconstruction_stencil_downstream_size )][k + z_varying * ( m - reconstruction_stencil_downstream_size )] * roe_eigenvector_left_temp[l];

                     // Compute characteristic advections to compute fluxes from left and right side of the face i+1/2
                     positive_characteristic_flux[m] += ( advection_characteristic_temp + alpha_fluxvectorsplitting[n] * u_characteristic_temp );
                     negative_characteristic_flux[m] += ( advection_characteristic_temp - alpha_fluxvectorsplitting[n] * u_characteristic_temp );
                  }
               }

               //apply WENO scheme to compute characteristic fluxes
               characteristic_flux[n] = 0.5 * ( SU::Reconstruction<ReconstructionStencil, SP::UpwindLeft>( positive_characteristic_flux, cell_size ) + SU::Reconstruction<ReconstructionStencil, SP::UpwindRight>( negative_characteristic_flux, cell_size ) );
            }

            // reconstruct the fluxes at face i+1/2
            auto const physical_flux = eigendecomposition_calculator_.TransformToPhysicalSpace( characteristic_flux, roe_eigenvectors_right[i_index][j_index][k_index] );

            for( unsigned int l = 0; l < MF::ANOE(); ++l ) {
               fluxes[l][i_index][j_index][k_index] += physical_flux[l];
            }
         }
      }
   }
}
