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
#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H
#include <array>

#include "block_definitions/block.h"
#include "solvers/eigendecomposition.h"
#include "materials/material_manager.h"

/**
 * @brief Interface to solve the underlying system of equations. Uses spatial reconstruction stencils to approximate the solution.
 */
template<typename DerivedRiemannSolver>
class RiemannSolver {

   friend DerivedRiemannSolver;

   EigenDecomposition const& eigendecomposition_calculator_;
   MaterialManager const& material_manager_;

   explicit RiemannSolver( MaterialManager const& material_manager, EigenDecomposition const& eigendecomposition_calculator ) : eigendecomposition_calculator_( eigendecomposition_calculator ),
                                                                                                                                material_manager_( material_manager ) {
      //Empty constructor besides initializer list.
   }

public:
   RiemannSolver()                       = delete;
   ~RiemannSolver()                      = default;
   RiemannSolver( RiemannSolver const& ) = delete;
   RiemannSolver& operator=( RiemannSolver const& ) = delete;
   RiemannSolver( RiemannSolver&& )                 = delete;
   RiemannSolver& operator=( RiemannSolver&& ) = delete;

   template<Direction DIR>
   std::array<double, MF::ANOE()> SolveRiemannProblem( MaterialName const material,
                                                       std::array<double, MF::ANOE()> const& state_face_left,
                                                       std::array<double, MF::ANOE()> const& state_face_right,
                                                       std::array<double, MF::ANOP()> const& prime_state_left,
                                                       std::array<double, MF::ANOP()> const& prime_state_right ) const {
      return static_cast<DerivedRiemannSolver const&>( *this ).template SolveRiemannProblemImplementation<DIR>( material, state_face_left, state_face_right, prime_state_left, prime_state_right );
   }
};

#endif// RIEMANN_SOLVER_H
