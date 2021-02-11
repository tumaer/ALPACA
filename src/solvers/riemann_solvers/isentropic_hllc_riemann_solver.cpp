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

#include "solvers/riemann_solvers/isentropic_hllc_riemann_solver.h"

#include "utilities/mathematical_functions.h"
#include "solvers/riemann_solvers/hll_signal_speed_calculator.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief Standard constructor using an already existing MaterialManager and EigenDecomposition object. See base class.
 * @param material_manager .
 * @param eigendecomposition_calculator .
 */
IsentropicHllcRiemannSolver::IsentropicHllcRiemannSolver( MaterialManager const& material_manager, EigenDecomposition const& eigendecomposition_calculator ) : RiemannSolver( material_manager, eigendecomposition_calculator ) {
   /* Empty besides initializer list*/
}

/**
 * @brief Solving the right hand side of the underlying system of equations. Using spatial reconstruction of cell averaged values (finite volume approach) and flux determination by HLLC procedure.
 * See base class.
 */
void IsentropicHllcRiemannSolver::UpdateImplementation( std::pair<MaterialName const, Block> const& mat_block, double const cell_size,
                                                        double ( &fluxes_x )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                                                        double ( &fluxes_y )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                                                        double ( &fluxes_z )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1] ) const {

   ComputeFluxes<Direction::X>( mat_block, fluxes_x, cell_size );

   if constexpr( CC::DIM() != Dimension::One ) {
      ComputeFluxes<Direction::Y>( mat_block, fluxes_y, cell_size );
   }

   if constexpr( CC::DIM() == Dimension::Three ) {
      ComputeFluxes<Direction::Z>( mat_block, fluxes_z, cell_size );
   }
}

/**
 * @brief Computes the cell face fluxes with the set stencil & HLLC procedure
 * @param mat_block The block and material information of the phase under consideration.
 * @param fluxes Reference to an array which is filled with the computed fluxes (indirect return parameter).
 * @param cell_size .
 * @tparam DIR Indicates which spatial direction is to be computed.
 * @note Hotpath function.
 */
template<Direction DIR>
void IsentropicHllcRiemannSolver::ComputeFluxes( std::pair<MaterialName const, Block> const& mat_block,
                                                 double ( &fluxes )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                                                 double const cell_size ) const {

   using ReconstructionStencil = ReconstructionStencilSetup::Concretize<reconstruction_stencil>::type;

   // declaration of applied arrays
   std::array<double, ReconstructionStencil::StencilSize()> reconstruction_array;// storage for reconstruction values
   std::array<double, MF::ANOE()> state_face_left;                               //variable vector containing interpolated states of left patch of cell face i/j/k+1/2
   std::array<double, MF::ANOE()> state_face_right;                              //variable vector containing interpolated states of right patch of cell face i/j/k+1/2
   std::array<double, MF::ANOE()> q_star_left;                                   // intermediate state vector between left-going wave and contact-interface
   std::array<double, MF::ANOE()> q_star_right;                                  // intermediate state vector between contact-interface and right-going wave
   std::array<double, MF::ANOE()> flux_left;                                     // F(state_face_left)
   std::array<double, MF::ANOE()> flux_right;                                    // F(state_face_right)

   constexpr unsigned int principal_momentum_index = ETI( MF::AME()[DTI( DIR )] );

   constexpr unsigned int x_start = DIR == Direction::X ? CC::FICX() - 1 : CC::FICX();
   constexpr unsigned int y_start = DIR == Direction::Y ? CC::FICY() - 1 : CC::FICY();
   constexpr unsigned int z_start = DIR == Direction::Z ? CC::FICZ() - 1 : CC::FICZ();

   constexpr unsigned int x_reconstruction_offset = DIR == Direction::X ? 1 : 0;
   constexpr unsigned int y_reconstruction_offset = DIR == Direction::Y ? 1 : 0;
   constexpr unsigned int z_reconstruction_offset = DIR == Direction::Z ? 1 : 0;

   constexpr unsigned int x_end = CC::LICX();
   constexpr unsigned int y_end = CC::LICY();
   constexpr unsigned int z_end = CC::LICZ();

   constexpr int total_to_internal_offset_x = CC::FICX() - 1;
   constexpr int total_to_internal_offset_y = CC::DIM() != Dimension::One ? static_cast<int>( CC::FICY() ) - 1 : -1;
   constexpr int total_to_internal_offset_z = CC::DIM() == Dimension::Three ? static_cast<int>( CC::FICZ() ) - 1 : -1;

   // Access the pair's elements directly.
   auto const& [material, block] = mat_block;

   // Required for Toro signal speeds
   double const gamma = material_manager_.GetMaterial( material ).GetEquationOfState().Gamma();

   for( unsigned int i = x_start; i <= x_end; ++i ) {
      for( unsigned int j = y_start; j <= y_end; ++j ) {
         for( unsigned int k = z_start; k <= z_end; ++k ) {
            // Shifted indices to match block index system and roe-ev index system
            int const i_index = i - total_to_internal_offset_x;
            int const j_index = j - total_to_internal_offset_y;
            int const k_index = k - total_to_internal_offset_z;

            // Reconstruct fluxes at face i+1/2 using a selected reconstruction procedure
            for( unsigned int n = 0; n < MF::ANOE(); ++n ) {
               for( unsigned int m = 0; m < ReconstructionStencil::StencilSize(); ++m ) {
                  reconstruction_array[m] = block.GetAverageBuffer( MF::ASOE()[n] )[i + x_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )]
                                                                                   [j + y_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )]
                                                                                   [k + z_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )];
               }
               state_face_left[n]  = SU::Reconstruction<ReconstructionStencil, SP::UpwindLeft>( reconstruction_array, cell_size );
               state_face_right[n] = SU::Reconstruction<ReconstructionStencil, SP::UpwindRight>( reconstruction_array, cell_size );
            }

            // Check for invalid cells due to ghost fluid method
            if( state_face_left[ETI( Equation::Mass )] <= std::numeric_limits<double>::epsilon() || state_face_right[ETI( Equation::Mass )] <= std::numeric_limits<double>::epsilon() ) continue;

            // Compute pressure, velocity and speed of sound for both cells for reconstructed values
            double const pressure_left  = material_manager_.GetMaterial( material ).GetEquationOfState().Pressure( state_face_left[ETI( Equation::Mass )], 0.0, 0.0, 0.0, 0.0 );
            double const pressure_right = material_manager_.GetMaterial( material ).GetEquationOfState().Pressure( state_face_right[ETI( Equation::Mass )], 0.0, 0.0, 0.0, 0.0 );

            // Check for invalid cells due to ghost fluid method
            if( pressure_left <= 0 || pressure_right <= 0 ) continue;

            double const velocity_left        = state_face_left[principal_momentum_index] / state_face_left[ETI( Equation::Mass )];
            double const velocity_right       = state_face_right[principal_momentum_index] / state_face_right[ETI( Equation::Mass )];
            double const speed_of_sound_left  = material_manager_.GetMaterial( material ).GetEquationOfState().SpeedOfSound( state_face_left[ETI( Equation::Mass )], pressure_left );
            double const speed_of_sound_right = material_manager_.GetMaterial( material ).GetEquationOfState().SpeedOfSound( state_face_right[ETI( Equation::Mass )], pressure_right );

            // Calculation of signal speeds
            auto const [wave_speed_left_simple, wave_speed_right_simple] = CalculateSignalSpeed( state_face_left[ETI( Equation::Mass )], state_face_right[ETI( Equation::Mass )],
                                                                                                 velocity_left, velocity_right,
                                                                                                 pressure_left, pressure_right,
                                                                                                 speed_of_sound_left, speed_of_sound_right,
                                                                                                 gamma );

            double const wave_speed_contact = ( ( pressure_right - pressure_left ) + ( state_face_left[ETI( Equation::Mass )] * velocity_left * ( wave_speed_left_simple - velocity_left ) - state_face_right[ETI( Equation::Mass )] * velocity_right * ( wave_speed_right_simple - velocity_right ) ) ) / ( state_face_left[ETI( Equation::Mass )] * ( wave_speed_left_simple - velocity_left ) - state_face_right[ETI( Equation::Mass )] * ( wave_speed_right_simple - velocity_right ) );
            double const wave_speed_left    = std::min( wave_speed_left_simple, 0.0 );
            double const wave_speed_right   = std::max( wave_speed_right_simple, 0.0 );

            double const chi_star_left  = ( wave_speed_left_simple - velocity_left ) / ( wave_speed_left_simple - wave_speed_contact );
            double const chi_star_right = ( wave_speed_right_simple - velocity_right ) / ( wave_speed_right_simple - wave_speed_contact );

            // Compute intermediate states (Toro 10.71 10.72 10.73) and calculate F(q)
            q_star_left[ETI( Equation::Mass )]     = state_face_left[ETI( Equation::Mass )] * chi_star_left;
            q_star_left[principal_momentum_index]  = state_face_left[ETI( Equation::Mass )] * chi_star_left * wave_speed_contact;
            q_star_right[ETI( Equation::Mass )]    = state_face_right[ETI( Equation::Mass )] * chi_star_right;
            q_star_right[principal_momentum_index] = state_face_right[ETI( Equation::Mass )] * chi_star_right * wave_speed_contact;

            flux_left[ETI( Equation::Mass )]     = state_face_left[principal_momentum_index];
            flux_left[principal_momentum_index]  = ( ( state_face_left[principal_momentum_index] * state_face_left[principal_momentum_index] ) / state_face_left[ETI( Equation::Mass )] ) + pressure_left;
            flux_right[ETI( Equation::Mass )]    = state_face_right[principal_momentum_index];
            flux_right[principal_momentum_index] = ( ( state_face_right[principal_momentum_index] * state_face_right[principal_momentum_index] ) / state_face_right[ETI( Equation::Mass )] ) + pressure_right;

            // minor momenta
            for( unsigned int d = 0; d < DTI( CC::DIM() ) - 1; ++d ) {
               // get the index of this minor momentum
               unsigned int const minor_momentum_index = ETI( MF::AME()[DTI( GetMinorDirection<DIR>( d ) )] );

               q_star_left[minor_momentum_index]  = state_face_left[ETI( Equation::Mass )] * chi_star_left * ( state_face_left[minor_momentum_index] / state_face_left[ETI( Equation::Mass )] );
               q_star_right[minor_momentum_index] = state_face_right[ETI( Equation::Mass )] * chi_star_right * ( state_face_right[minor_momentum_index] / state_face_right[ETI( Equation::Mass )] );

               flux_left[minor_momentum_index]  = velocity_left * state_face_left[minor_momentum_index];
               flux_right[minor_momentum_index] = velocity_right * state_face_right[minor_momentum_index];
            }

            for( unsigned int n = 0; n < MF::ANOE(); ++n ) {
               double const flux_star_left  = flux_left[n] + ( wave_speed_left * ( q_star_left[n] - state_face_left[n] ) );
               double const flux_star_right = flux_right[n] + ( wave_speed_right * ( q_star_right[n] - state_face_right[n] ) );
               // Attention plus used for fluxes as other terms will also increment elements.
               fluxes[n][i_index][j_index][k_index] += ( 0.5 * ( 1.0 + Signum( wave_speed_contact ) ) * flux_star_left + 0.5 * ( 1.0 - Signum( wave_speed_contact ) ) * flux_star_right );
            }
         }// k
      }   // j
   }      // i
}
