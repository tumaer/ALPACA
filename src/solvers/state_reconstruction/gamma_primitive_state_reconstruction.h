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
#ifndef GAMMA_PRIMITIVE_STATE_RECONSTRUCTION_H
#define GAMMA_PRIMITIVE_STATE_RECONSTRUCTION_H

#include "solvers/state_reconstruction/state_reconstruction.h"
#include "prime_states/prime_state_handler.h"
#include "block_definitions/field_material_definitions.h"

#include "stencils/stencil_utilities.h"

/**
 * @brief Discretization of the spatial reconstruction scheme using primitive states for reconstruction ( except for gamma and pi where always the conservative quantity is used ).
 */
class GammaPrimitiveStateReconstruction : public StateReconstruction<GammaPrimitiveStateReconstruction> {

   friend StateReconstruction;

   /**
    * @brief Procedure to reconstruct the conservatives/primitive states at cell faces.
    * @tparam DIR spatial direction the reconstruction has to be performed.
    * @param block Block of the phase under consideration.
    * @param eos Underlying equation of state of the phase under consideration used to convert primes and conservatives.
    * @param roe_eigenvectors_left .
    * @param roe_eigenvectors_right .
    * @param cell_size .
    * @param i .
    * @param k .
    * @param j .
    * @return tuple containing left and right reconstructed primitive and conservative states.
    */
   template<Direction DIR, ReconstructionStencils RECON>
   std::tuple<std::array<double, MF::ANOE()>, std::array<double, MF::ANOE()>, std::array<double, MF::ANOP()>, std::array<double, MF::ANOP()>>
   SolveStateReconstructionImplementation( Block const& block,
                                           EquationOfState const& eos,
                                           double const ( & )[MF::ANOE()][MF::ANOE()],
                                           double const ( & )[MF::ANOE()][MF::ANOE()],
                                           double const cell_size,
                                           unsigned int const i,
                                           unsigned int const j,
                                           unsigned int const k ) const {

      using ReconstructionStencil                    = typename ReconstructionStencilSetup::Concretize<RECON>::type;
      constexpr unsigned int x_reconstruction_offset = DIR == Direction::X ? 1 : 0;
      constexpr unsigned int y_reconstruction_offset = DIR == Direction::Y ? 1 : 0;
      constexpr unsigned int z_reconstruction_offset = DIR == Direction::Z ? 1 : 0;
      std::array<double, MF::ANOP()> reconstructed_primes_minus;
      std::array<double, MF::ANOP()> reconstructed_primes_plus;

      std::array<double, ReconstructionStencil::StencilSize()> reconstruction_array;
      std::array<double, MF::ANOE()> reconstructed_conservatives_minus;
      std::array<double, MF::ANOE()> reconstructed_conservatives_plus;
      for( unsigned int n = 0; n < MF::ANOP(); ++n ) {
         for( unsigned int m = 0; m < ReconstructionStencil::StencilSize(); ++m ) {
            reconstruction_array[m] = block.GetPrimeStateBuffer( MF::ASOP()[n] )[i + x_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )]
                                                                                [j + y_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )]
                                                                                [k + z_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )];
         }// M-Loop

         reconstructed_primes_minus[n] = SU::Reconstruction<ReconstructionStencil, SP::UpwindLeft>( reconstruction_array, cell_size );
         reconstructed_primes_plus[n]  = SU::Reconstruction<ReconstructionStencil, SP::UpwindRight>( reconstruction_array, cell_size );
      }// N-Loop

      for( unsigned int m = 0; m < ReconstructionStencil::StencilSize(); ++m ) {
         reconstruction_array[m] = block.GetAverageBuffer( MF::ASOE()[ETI( Equation::Gamma )] )[i + x_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )]
                                                                                               [j + y_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )]
                                                                                               [k + z_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )];
      }
      reconstructed_conservatives_minus[ETI( Equation::Gamma )] = SU::Reconstruction<ReconstructionStencil, SP::UpwindLeft>( reconstruction_array, cell_size );
      reconstructed_conservatives_plus[ETI( Equation::Gamma )]  = SU::Reconstruction<ReconstructionStencil, SP::UpwindRight>( reconstruction_array, cell_size );

      for( unsigned int m = 0; m < ReconstructionStencil::StencilSize(); ++m ) {
         reconstruction_array[m] = block.GetAverageBuffer( MF::ASOE()[ETI( Equation::Pi )] )[i + x_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )]
                                                                                            [j + y_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )]
                                                                                            [k + z_reconstruction_offset * ( m - ReconstructionStencil::DownstreamStencilSize() )];
      }
      reconstructed_conservatives_minus[ETI( Equation::Pi )] = SU::Reconstruction<ReconstructionStencil, SP::UpwindLeft>( reconstruction_array, cell_size );
      reconstructed_conservatives_plus[ETI( Equation::Pi )]  = SU::Reconstruction<ReconstructionStencil, SP::UpwindRight>( reconstruction_array, cell_size );

      reconstructed_primes_minus[PTI( PrimeState::gamma )] = GammaModelStiffenedGas::CalculatePrimeGamma( reconstructed_conservatives_minus[ETI( Equation::Gamma )] );
      reconstructed_primes_plus[PTI( PrimeState::gamma )]  = GammaModelStiffenedGas::CalculatePrimeGamma( reconstructed_conservatives_plus[ETI( Equation::Gamma )] );
      reconstructed_primes_minus[PTI( PrimeState::pi )]    = GammaModelStiffenedGas::CalculatePrimePi( reconstructed_primes_minus[PTI( PrimeState::gamma )], reconstructed_conservatives_minus[ETI( Equation::Pi )] );
      reconstructed_primes_plus[PTI( PrimeState::pi )]     = GammaModelStiffenedGas::CalculatePrimePi( reconstructed_primes_plus[PTI( PrimeState::gamma )], reconstructed_conservatives_plus[ETI( Equation::Pi )] );

      PrimeStatesToConservatives( eos, reconstructed_primes_minus, reconstructed_conservatives_minus );
      PrimeStatesToConservatives( eos, reconstructed_primes_plus, reconstructed_conservatives_plus );
      return std::make_tuple( reconstructed_conservatives_minus, reconstructed_conservatives_plus, reconstructed_primes_minus, reconstructed_primes_plus );
   }

public:
   GammaPrimitiveStateReconstruction() : StateReconstruction() {}
   ~GammaPrimitiveStateReconstruction()                                          = default;
   GammaPrimitiveStateReconstruction( GammaPrimitiveStateReconstruction const& ) = delete;
   GammaPrimitiveStateReconstruction& operator=( GammaPrimitiveStateReconstruction const& ) = delete;
   GammaPrimitiveStateReconstruction( GammaPrimitiveStateReconstruction&& )                 = delete;
   GammaPrimitiveStateReconstruction& operator=( GammaPrimitiveStateReconstruction&& ) = delete;
};

#endif// GAMMA_PRIMITIVE_STATE_RECONSTRUCTION_H
