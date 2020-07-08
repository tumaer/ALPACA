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
#include "exact_iterative_interface_riemann_solver.h"

/**
 * @brief      Calculates the shock/rarefaction relation and its derivative. See \cite Toro2009.
 *
 * @param[in]  initial_root           The interface pressure of the current iteration.
 * @param[in]  p                      The left/right pressure defining the interface Riemann problem.
 * @param[in]  pressure_function      A constant dependent on the EoS.
 * @param[in]  one_pressure_function  The inverse of the pressure function.
 * @param[in]  pressure_constant      The pressure constant B.
 * @param[in]  A                      A constant. See IterationConstants::A() for details.
 * @param[in]  B                      A constant. See IterationConstants::B() for details.
 * @param[in]  C                      A constant. See IterationConstants::C() for details.
 * @param[in]  D                      A constant. See IterationConstants::D() for details.
 *
 * @return     The shock/rarefaction relation and its derivative.
 */
std::array<double, 2> ExactIterativeInterfaceRiemannSolver::ObtainFunctionAndDerivativeImplementation( double const initial_root, double const p,
                                                                                                       double const pressure_function, double const one_pressure_function, double const pressure_constant, double const A, double const B, double const C, double const D ) const {

   if( initial_root > p ) {
      return { ShockRelations::Function( IterationUtilities::MaterialPressureFunction( initial_root, pressure_constant ), pressure_function, A, B ),
               ShockRelations::Derivative( IterationUtilities::MaterialPressureFunction( initial_root, pressure_constant ), pressure_function, A, B ) };
   } else {
      return { RarefactionRelations::Function( IterationUtilities::MaterialPressureFunction( initial_root, pressure_constant ), one_pressure_function, C, D ),
               RarefactionRelations::Derivative( IterationUtilities::MaterialPressureFunction( initial_root, pressure_constant ), one_pressure_function, C, D ) };
   }
}