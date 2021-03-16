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
#include "solvers/convective_term_contributions/finite_volume_scheme.h"
#include "stencils/stencil_utilities.h"
#include "solvers/state_reconstruction.h"

static_assert( !( active_equations == EquationSet::Isentropic && state_reconstruction_type == StateReconstructionType::RoeCharacteristic ), "RoeCharacteristic reconstruction not implemented for isentropic equations!" );

/**
 * @brief Standard constructor using an already existing MaterialManager and EigenDecomposition object.
 * @param material_manager .
 * @param eigendecomposition_calculator .
 */
FiniteVolumeScheme::FiniteVolumeScheme( MaterialManager const& material_manager, EigenDecomposition const& eigendecomposition_calculator ) : ConvectiveTermSolver( material_manager, eigendecomposition_calculator ),
                                                                                                                                             riemann_solver_( material_manager, eigendecomposition_calculator_ ) {
   /* Empty besides initializer list*/
}

/**
 * @brief Solving the convective term of the system. Using dimension splitting for fluxes in x, y, and z- direction. Also See base class.
 * @note Hotpath function.
 */
void FiniteVolumeScheme::UpdateImplementation(
      std::pair<MaterialName const, Block> const& mat_block, double const cell_size,
      double ( &fluxes_x )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double ( &fluxes_y )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double ( &fluxes_z )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1] ) const {

   double roe_eigenvectors_left[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()];
   double roe_eigenvectors_right[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()];
   double roe_eigenvalues[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()];

   constexpr bool require_eigendecomposition = ( active_equations == EquationSet::NavierStokes || active_equations == EquationSet::Euler ) && state_reconstruction_type == StateReconstructionType::RoeCharacteristic;

   if constexpr( require_eigendecomposition ) {
      // Setting the eigenvectors, eigenvalues to zero is necessary for two-phase simulations as not every entry is necessarily set during eigendecomposition
      for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
         for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
            for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
               for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
                  for( unsigned int f = 0; f < MF::ANOE(); ++f ) {
                     roe_eigenvectors_left[i][j][k][e][f]  = 0.0;
                     roe_eigenvectors_right[i][j][k][e][f] = 0.0;
                  }
                  roe_eigenvalues[i][j][k][e] = 0.0;
               }
            }
         }
      }

      eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::X>( mat_block, roe_eigenvectors_left, roe_eigenvectors_right, roe_eigenvalues );
   }
   ComputeFluxes<Direction::X>( mat_block, fluxes_x, roe_eigenvectors_left, roe_eigenvectors_right, cell_size );

   if constexpr( CC::DIM() != Dimension::One ) {
      if constexpr( require_eigendecomposition ) {
         eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::Y>( mat_block, roe_eigenvectors_left, roe_eigenvectors_right, roe_eigenvalues );
      }
      ComputeFluxes<Direction::Y>( mat_block, fluxes_y, roe_eigenvectors_left, roe_eigenvectors_right, cell_size );
   }

   if constexpr( CC::DIM() == Dimension::Three ) {
      if constexpr( require_eigendecomposition ) {
         eigendecomposition_calculator_.ComputeRoeEigendecomposition<Direction::Z>( mat_block, roe_eigenvectors_left, roe_eigenvectors_right, roe_eigenvalues );
      }
      ComputeFluxes<Direction::Z>( mat_block, fluxes_z, roe_eigenvectors_left, roe_eigenvectors_right, cell_size );
   }
}

/**
 * @brief Computes the convective cell face fluxes by performing a finite-volume state reconstruction and solving a Riemann problem afterwards.
 * @param mat_block The block and material information of the phase under consideration.
 * @param fluxes Reference to an array which is filled with the computed fluxes (indirect return parameter).
 * @param roe_eigenvectors_left .
 * @param roe_eigenvectors_right .
 * @param cell_size .
 * @tparam DIR Indicates which spatial direction is to be computed.
 * @note Hotpath function.
 */
template<Direction DIR>
void FiniteVolumeScheme::ComputeFluxes( std::pair<MaterialName const, Block> const& mat_block,
                                        double ( &fluxes )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                                        double const ( &Roe_eigenvectors_left )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()],
                                        double const ( &Roe_eigenvectors_right )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()],
                                        double const cell_size ) const {

   constexpr unsigned int x_start = DIR == Direction::X ? CC::FICX() - 1 : CC::FICX();
   constexpr unsigned int y_start = DIR == Direction::Y ? CC::FICY() - 1 : CC::FICY();
   constexpr unsigned int z_start = DIR == Direction::Z ? CC::FICZ() - 1 : CC::FICZ();

   constexpr unsigned int x_end = CC::LICX();
   constexpr unsigned int y_end = CC::LICY();
   constexpr unsigned int z_end = CC::LICZ();

   constexpr int total_to_internal_offset_x = CC::FICX() - 1;
   constexpr int total_to_internal_offset_y = CC::DIM() != Dimension::One ? static_cast<int>( CC::FICY() ) - 1 : -1;
   constexpr int total_to_internal_offset_z = CC::DIM() == Dimension::Three ? static_cast<int>( CC::FICZ() ) - 1 : -1;

   // Access the pair's elements directly.
   auto const& [material, block] = mat_block;

   for( unsigned int i = x_start; i <= x_end; ++i ) {
      for( unsigned int j = y_start; j <= y_end; ++j ) {
         for( unsigned int k = z_start; k <= z_end; ++k ) {
            // Shifted indices to match block index system and roe-ev index system
            int const i_index = i - total_to_internal_offset_x;
            int const j_index = j - total_to_internal_offset_y;
            int const k_index = k - total_to_internal_offset_z;

            auto const [reconstructed_conservatives_left, reconstructed_conservatives_right,
                        reconstructed_primes_left, reconstructed_primes_right] = StateReconstruction<DIR, reconstruction_stencil>( block,
                                                                                                                                   material_manager_.GetMaterial( material ).GetEquationOfState(),
                                                                                                                                   Roe_eigenvectors_left[i_index][j_index][k_index],
                                                                                                                                   Roe_eigenvectors_right[i_index][j_index][k_index],
                                                                                                                                   cell_size, i, j, k );
            // To check for invalid cells due to ghost fluid method
            double const B = active_equations == EquationSet::Isentropic ? 0.0 : material_manager_.GetMaterial( material ).GetEquationOfState().B();
            if( reconstructed_conservatives_left[ETI( Equation::Mass )] <= std::numeric_limits<double>::epsilon() || reconstructed_conservatives_right[ETI( Equation::Mass )] <= std::numeric_limits<double>::epsilon() ) continue;
            if( reconstructed_primes_left[PTI( PrimeState::Pressure )] <= -B || reconstructed_primes_right[PTI( PrimeState::Pressure )] <= -B ) continue;

            auto const face_fluxes = riemann_solver_.SolveRiemannProblem<DIR>( material, reconstructed_conservatives_left, reconstructed_conservatives_right, reconstructed_primes_left, reconstructed_primes_right );

            for( unsigned int n = 0; n < MF::ANOE(); ++n ) {
               fluxes[n][i_index][j_index][k_index] += face_fluxes[n];
            }
         }// k
      }   // j
   }      // i
}
