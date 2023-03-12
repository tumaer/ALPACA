//===----------------------- numerical_setup.h ----------------------------===//
//
//                                 ALPACA
//
// Part of ALPACA, under the GNU General Public License as published by
// the Free Software Foundation version 3.
// SPDX-License-Identifier: GPL-3.0-only
//
// If using this code in an academic setting, please cite the following:
// @article{hoppe2022parallel,
//  title={A parallel modular computing environment for three-dimensional
//  multiresolution simulations of compressible flows},
//  author={Hoppe, Nils and Adami, Stefan and Adams, Nikolaus A},
//  journal={Computer Methods in Applied Mechanics and Engineering},
//  volume={391},
//  pages={114486},
//  year={2022},
//  publisher={Elsevier}
// }
//
//===----------------------------------------------------------------------===//
#ifndef NUMERICAL_SETUP_H
#define NUMERICAL_SETUP_H

// TIME_INTEGRATION_SCHEME
enum class TimeIntegrators { RK2, RK3 };
constexpr TimeIntegrators time_integrator = TimeIntegrators::RK3;

// MULTI_PHASE_MANAGER
enum class PhaseManagers { TwoPhase, NumerousPhase };
constexpr PhaseManagers phase_manager = PhaseManagers::TwoPhase;

// LEVELSET_ADVECTOR
enum class LevelsetAdvectors {
  DerivativeStencil,
  ReconstructionStencil,
  HjReconstructionStencil,
  HjDerivativeStencil
};
constexpr LevelsetAdvectors levelset_advector =
    LevelsetAdvectors::HjReconstructionStencil;

// LEVELSET_REINITIALIZER
enum class LevelsetReinitializers { Min, Weno, Explicit };
constexpr LevelsetReinitializers levelset_reinitializer =
    LevelsetReinitializers::Weno;

// MIXING_METHOD
enum class CutCellMixers { ApertureBased, NormalBased, Lauer };
constexpr CutCellMixers cut_cell_mixer = CutCellMixers::ApertureBased;

// EXTENSION_METHOD
enum class Extenders { Fedkiw, Upwind, Explicit };
constexpr Extenders extender = Extenders::Fedkiw;

// GEOMETRY_CALCULATOR
enum class GeometryCalculators { MarchingCubes };
constexpr GeometryCalculators geometry_calculator =
    GeometryCalculators::MarchingCubes;

// INTERFACE_RIEMANN_SOLVER
enum class InterfaceRiemannSolvers { Linearized, Exact, TwoRarefaction, Hllc };
constexpr InterfaceRiemannSolvers interface_riemann_solver =
    InterfaceRiemannSolvers::Linearized;

// SCALE_SEPARATOR
enum class ScaleSeparators { TwoPhase };
constexpr ScaleSeparators scale_separator = ScaleSeparators::TwoPhase;

// INTERFACE_EXTENDER
enum class InterfaceExtenders { TwoPhase };
constexpr InterfaceExtenders interface_extender = InterfaceExtenders::TwoPhase;

// BUFFER_HANDLER
enum class BufferHandlers { TwoPhase };
constexpr BufferHandlers buffer_handler = BufferHandlers::TwoPhase;

#endif // NUMERICAL_SETUP_H
