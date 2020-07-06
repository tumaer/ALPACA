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
#ifndef SPACE_SOLVER_H
#define SPACE_SOLVER_H

#include "user_specifications/numerical_setup.h"
#include "block_definitions/block.h"
#include "topology/node.h"
#include "materials/material_manager.h"
#include "source_term_solver.h"
#include "eigendecomposition.h"
#include "interface_interaction/interface_term_solver.h"

#include "levelset/levelset_advector/levelset_advector_setup.h"
#include "solvers/riemann_solvers/riemann_solver_setup.h"

using RiemannSolverConcretization = RiemannSolverSetup::Concretize<riemann_solver>::type;
using LevelsetAdvectorConcretization = LevelsetAdvectorSetup::Concretize<levelset_advector>::type;

/**
 * @brief The SpaceSolver solves right side of the underlying system of equations (including source terms) using a Riemann solver of choice with an interchangable stencil.
 */
class SpaceSolver {

   EigenDecomposition const eigendecomposition_calculator_;
   RiemannSolverConcretization const riemann_solver_;
   SourceTermSolver const source_term_solver_;
   InterfaceTermSolver const interface_term_solver_;
   LevelsetAdvectorConcretization const levelset_advector_;

public:
   SpaceSolver() = delete;
   explicit SpaceSolver( MaterialManager const& material_manager, std::array<double, 3> gravity);
   ~SpaceSolver() = default;
   SpaceSolver( SpaceSolver const& ) = delete;
   SpaceSolver& operator=( SpaceSolver const& ) = delete;
   SpaceSolver( SpaceSolver&& ) = delete;
   SpaceSolver& operator=( SpaceSolver&& ) = delete;

   void UpdateFluxes( Node& node ) const;
   void ComputeMaxEigenvaluesForPhase( std::pair<MaterialName const, Block> const& mat_block, double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()] ) const;
   void SetFluxFunctionGlobalEigenvalues( double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()] ) const;
};

#endif // SPACE_SOLVER_H
