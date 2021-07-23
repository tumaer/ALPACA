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
#ifndef RIEMANN_SOLVER_SETTINGS_H
#define RIEMANN_SOLVER_SETTINGS_H

#include "enums/flux_splitting.h"
#include "enums/signal_speed.h"

enum class ConvectiveTermSolvers { FluxSplitting,
                                   FiniteVolume };
constexpr ConvectiveTermSolvers convective_term_solver = ConvectiveTermSolvers::FluxSplitting;

namespace FluxSplittingSettings {
   /* FluxSplitting options are:
    * Roe | LocalLaxFriedrichs | GlobalLaxFriedrichs | Roe_M | LocalLaxFriedrichs_M
    * Roe_M and LocalLaxFriedrichs_M according to \cite Fleischmann20
    */
   constexpr FluxSplitting flux_splitting_scheme = FluxSplitting::Roe;

   /* Phi in \cite Fleischmann20.
    * Limits the speed of sound in the eigenvalue calculation of Roe-M and LLF-M
    */
   constexpr double low_mach_number_limit_factor = 5.0;
}// namespace FluxSplittingSettings

namespace FiniteVolumeSettings {
   // RIEMANN_SOLVER
   enum class RiemannSolvers { Hllc,
                               Hllc_LM,
                               Hll };
   constexpr RiemannSolvers riemann_solver = RiemannSolvers::Hllc;

   /* Signal Speed choices for HLL-type solvers are:
    * Einfeldt \cite Einfeldt88 | Davis \cite Davis88 | Toro \cite Toro94 | Arithmetic \cite Coralic14
    */
   constexpr SignalSpeed signal_speed_selection = SignalSpeed::Einfeldt;
}// namespace FiniteVolumeSettings

#endif// RIEMANN_SOLVER_SETTINGS_H
