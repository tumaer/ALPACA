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
#ifndef NUMERICAL_SETUP_H
#define NUMERICAL_SETUP_H

// TIME_INTEGRATION_SCHEME
enum class TimeIntegrators { RK2,
                             RK3 };
constexpr TimeIntegrators time_integrator = TimeIntegrators::RK3;

// MULTI_PHASE_MANAGER
enum class PhaseManagers { TwoPhase,
                           NumerousPhase };
constexpr PhaseManagers phase_manager = PhaseManagers::TwoPhase;

// LEVELSET_ADVECTOR
enum class LevelsetAdvectors { DerivativeStencil,
                               ReconstructionStencil,
                               HjReconstructionStencil,
                               HjDerivativeStencil };
constexpr LevelsetAdvectors levelset_advector = LevelsetAdvectors::HjReconstructionStencil;

// LEVELSET_REINITIALIZER
enum class LevelsetReinitializers { Min,
                                    Weno,
                                    Explicit };
constexpr LevelsetReinitializers levelset_reinitializer = LevelsetReinitializers::Weno;

// MIXING_METHOD
enum class CutCellMixers { ApertureBased,
                           NormalBased,
                           Lauer };
constexpr CutCellMixers cut_cell_mixer = CutCellMixers::ApertureBased;

// EXTENSION_METHOD
enum class Extenders { Fedkiw,
                       Upwind,
                       Explicit };
constexpr Extenders extender = Extenders::Fedkiw;

// GEOMETRY_CALCULATOR
enum class GeometryCalculators { MarchingCubes };
constexpr GeometryCalculators geometry_calculator = GeometryCalculators::MarchingCubes;

// INTERFACE_RIEMANN_SOLVER
enum class InterfaceRiemannSolvers { Linearized,
                                     Exact,
                                     TwoRarefaction,
                                     Hllc };
constexpr InterfaceRiemannSolvers interface_riemann_solver = InterfaceRiemannSolvers::Linearized;

// SCALE_SEPARATOR
enum class ScaleSeparators { TwoPhase };
constexpr ScaleSeparators scale_separator = ScaleSeparators::TwoPhase;

// INTERFACE_EXTENDER
enum class InterfaceExtenders { TwoPhase };
constexpr InterfaceExtenders interface_extender = InterfaceExtenders::TwoPhase;

// BUFFER_HANDLER
enum class BufferHandlers { TwoPhase };
constexpr BufferHandlers buffer_handler = BufferHandlers::TwoPhase;

#endif// NUMERICAL_SETUP_H
