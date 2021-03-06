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
#ifndef FLUX_SPLITTING_SCHEME_H
#define FLUX_SPLITTING_SCHEME_H

#include "solvers/convective_term_contributions/convective_term_solver.h"
#include "block_definitions/block.h"
#include "enums/direction_definition.h"
#include "utilities/helper_functions.h"
#include "materials/equation_of_state.h"
#include "materials/material_manager.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Discretization of the convective term solver using a flux-splitting procedure.
 */
class FluxSplittingScheme : public ConvectiveTermSolver<FluxSplittingScheme> {

   friend ConvectiveTermSolver;

   template<Direction DIR>
   void ComputeFluxes( Block const& b, double ( &fluxes )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                       double ( &advection )[MF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()], double const cell_size,
                       double ( &roe_eigenvectors_left )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()],
                       double ( &roe_eigenvectors_right )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()][MF::ANOE()],
                       double ( &fluxfunction_wavespeed )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][MF::ANOE()] ) const;

   void UpdateImplementation( std::pair<MaterialName const, Block> const& mat_block, double const cell_size,
                              double ( &fluxes_x )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                              double ( &fluxes_y )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                              double ( &fluxes_z )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                              double ( &volume_forces )[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()] ) const;

public:
   FluxSplittingScheme() = delete;
   explicit FluxSplittingScheme( MaterialManager const& material_manager, EigenDecomposition const& eigendecomposition_calculator );
   ~FluxSplittingScheme()                            = default;
   FluxSplittingScheme( FluxSplittingScheme const& ) = delete;
   FluxSplittingScheme& operator=( FluxSplittingScheme const& ) = delete;
   FluxSplittingScheme( FluxSplittingScheme&& )                 = delete;
   FluxSplittingScheme& operator=( FluxSplittingScheme&& ) = delete;
};

#endif// FLUX_SPLITTING_SCHEME_H
