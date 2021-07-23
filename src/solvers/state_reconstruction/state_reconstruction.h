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
#ifndef STATE_RECONSTRUCTION_H
#define STATE_RECONSTRUCTION_H

#include "user_specifications/stencil_setup.h"

/**
 * @brief Helper function to create the index sequence used to enforce symmetry while summing up conservative equation contributions in characteristic decomposition.
 * @tparam RemainingIndices Zero-based index sequence representing the non-momentum equations. Are transformed into real equation indices.
 * @return The created index sequence.
 */
template<std::size_t... NonMomentumIndices>
constexpr std::array<std::array<unsigned int, MF::ANOE()>, DTI( CC::DIM() )> MakeConservativeEquationSummationSequence( std::index_sequence<NonMomentumIndices...> const ) {
#if DIMENSION == 1
   return { { { ETI( Equation::MomentumX ), ETI( MF::NthNonMomentumEquation( NonMomentumIndices ) )... } } };
#elif DIMENSION == 2
   return { { { ETI( Equation::MomentumX ), ETI( Equation::MomentumY ), ETI( MF::NthNonMomentumEquation( NonMomentumIndices ) )... },
              { ETI( Equation::MomentumX ), ETI( Equation::MomentumY ), ETI( MF::NthNonMomentumEquation( NonMomentumIndices ) )... } } };
#else
   return { { { ETI( Equation::MomentumY ), ETI( Equation::MomentumZ ), ETI( Equation::MomentumX ), ETI( MF::NthNonMomentumEquation( NonMomentumIndices ) )... },
              { ETI( Equation::MomentumZ ), ETI( Equation::MomentumX ), ETI( Equation::MomentumY ), ETI( MF::NthNonMomentumEquation( NonMomentumIndices ) )... },
              { ETI( Equation::MomentumX ), ETI( Equation::MomentumY ), ETI( Equation::MomentumZ ), ETI( MF::NthNonMomentumEquation( NonMomentumIndices ) )... } } };
#endif
}

/**
 * @brief Interface to apply the spatial reconsturction scheme.
 */
template<typename DerivedStateReconstruction>
class StateReconstruction {

   friend DerivedStateReconstruction;

   explicit constexpr StateReconstruction() = default;

public:
   ~StateReconstruction()                            = default;
   StateReconstruction( StateReconstruction const& ) = delete;
   StateReconstruction& operator=( StateReconstruction const& ) = delete;
   StateReconstruction( StateReconstruction&& )                 = delete;
   StateReconstruction& operator=( StateReconstruction&& ) = delete;

   /**
    * @brief Solves the first-order Riemann problem using left and right state vectors.
    * @tparam direction.
    * @tparam reconstruction scheme.
    * @param block block.
    * @param eos applied equation of state.
    * @param Roe_eigenvectors_left, Roe_eigenvectors_right Left and right Roe eigenvector matrices for characteristic decomposition
    * @param cell_size cell size.
    * @param i,j,k cell indizes in block.
    * @return The reconstructed left/right conservative and primitive states.
    */
   template<Direction DIR, ReconstructionStencils RECON>
   std::tuple<std::array<double, MF::ANOE()>, std::array<double, MF::ANOE()>, std::array<double, MF::ANOP()>, std::array<double, MF::ANOP()>> SolveStateReconstruction( Block const& block,
                                                                                                                                                                        EquationOfState const& eos,
                                                                                                                                                                        double const ( &Roe_eigenvectors_left )[MF::ANOE()][MF::ANOE()],
                                                                                                                                                                        double const ( &Roe_eigenvectors_right )[MF::ANOE()][MF::ANOE()],
                                                                                                                                                                        double const cell_size,
                                                                                                                                                                        unsigned int const i,
                                                                                                                                                                        unsigned int const j,
                                                                                                                                                                        unsigned int const k ) const {
      return static_cast<DerivedStateReconstruction const&>( *this ).template SolveStateReconstructionImplementation<DIR, RECON>( block, eos, Roe_eigenvectors_left, Roe_eigenvectors_right, cell_size, i, j, k );
   }
};

#endif// STATE_RECONSTRUCTION_H
