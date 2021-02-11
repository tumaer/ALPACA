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

#ifndef ISENTROPIC_HLL_RIEMANN_SOLVER_H
#define ISENTROPIC_HLL_RIEMANN_SOLVER_H

#include "solvers/riemann_solvers/riemann_solver.h"
#include "block_definitions/block.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Discretization of the Riemann solver using the HLL procedure derived from \cite Toro2009, chapter 10.3 for isentropic cases.
 */
class IsentropicHllRiemannSolver : public RiemannSolver<IsentropicHllRiemannSolver> {

   friend RiemannSolver;

   // is used to distinguish principal and secondary momenta such that one single Riemann solver
   // routine can be used for all three spatial directions
   static constexpr std::array<std::array<unsigned int, 3>, 3> momentum_order_ = { { { ETI( Equation::MomentumX ), ETI( Equation::MomentumY ), ETI( Equation::MomentumZ ) },
                                                                                     { ETI( Equation::MomentumY ), ETI( Equation::MomentumX ), ETI( Equation::MomentumZ ) },
                                                                                     { ETI( Equation::MomentumZ ), ETI( Equation::MomentumX ), ETI( Equation::MomentumY ) } } };

   template<Direction DIR>
   void ComputeFluxes( std::pair<MaterialName const, Block> const& mat_block, double ( &fluxes )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                       double const cell_size ) const;

   void UpdateImplementation( std::pair<MaterialName const, Block> const& mat_block, double const cell_size,
                              double ( &fluxes_x )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                              double ( &fluxes_y )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                              double ( &fluxes_z )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1] ) const;

public:
   IsentropicHllRiemannSolver() = delete;
   explicit IsentropicHllRiemannSolver( MaterialManager const& material_manager, EigenDecomposition const& eigendecomposition_calculator );
   ~IsentropicHllRiemannSolver()                                   = default;
   IsentropicHllRiemannSolver( IsentropicHllRiemannSolver const& ) = delete;
   IsentropicHllRiemannSolver& operator=( IsentropicHllRiemannSolver const& ) = delete;
   IsentropicHllRiemannSolver( IsentropicHllRiemannSolver&& )                 = delete;
   IsentropicHllRiemannSolver& operator=( IsentropicHllRiemannSolver&& ) = delete;
};

#endif// HLL_RIEMANN_SOLVER_H
