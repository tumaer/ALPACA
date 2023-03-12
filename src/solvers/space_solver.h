//===------------------------- space_solver.h -----------------------------===//
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
#ifndef SPACE_SOLVER_H
#define SPACE_SOLVER_H

#include "block_definitions/block.h"
#include "eigendecomposition.h"
#include "interface_interaction/interface_term_solver.h"
#include "materials/material_manager.h"
#include "source_term_solver.h"
#include "topology/node.h"
#include "user_specifications/numerical_setup.h"

#include "levelset/levelset_advector/levelset_advector_setup.h"
#include "solvers/convective_term_contributions/convective_term_solver_setup.h"

using ConvectiveTermSolverConcretization =
    ConvectiveTermSolverSetup::Concretize<convective_term_solver>::type;
using LevelsetAdvectorConcretization =
    LevelsetAdvectorSetup::Concretize<levelset_advector>::type;

static_assert(!(active_equations == EquationSet::GammaModel &&
                CC::Axisymmetric()),
              "Axisymmetric terms are not implemented for Gamma-Model");

/**
 * @brief The SpaceSolver solves right side of the underlying system of
 * equations (including source terms) using a Riemann solver of choice with an
 * interchangable spatial stencil.
 */
class SpaceSolver {

  EigenDecomposition const eigendecomposition_calculator_;
  ConvectiveTermSolverConcretization const convective_term_solver_;
  SourceTermSolver const source_term_solver_;
  InterfaceTermSolver const interface_term_solver_;
  MaterialManager const &material_manager_;
  LevelsetAdvectorConcretization const levelset_advector_;

public:
  SpaceSolver() = delete;
  explicit SpaceSolver(MaterialManager const &material_manager,
                       std::array<double, 3> gravity);
  ~SpaceSolver() = default;
  SpaceSolver(SpaceSolver const &) = delete;
  SpaceSolver &operator=(SpaceSolver const &) = delete;
  SpaceSolver(SpaceSolver &&) = delete;
  SpaceSolver &operator=(SpaceSolver &&) = delete;

  void UpdateFluxes(Node &node) const;
  void UpdateLevelsetFluxes(Node &node) const;
  void ComputeMaxEigenvaluesForPhase(
      std::pair<MaterialName const, Block> const &mat_block,
      double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()]) const;
  void SetFluxFunctionGlobalEigenvalues(
      double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()]) const;
};

#endif // SPACE_SOLVER_H
