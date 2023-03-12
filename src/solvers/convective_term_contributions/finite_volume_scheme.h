//===----------------------- finite_volume_scheme.h -----------------------===//
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
#ifndef FINITE_VOLUME_SCHEME_H
#define FINITE_VOLUME_SCHEME_H

#include "block_definitions/block.h"
#include "enums/direction_definition.h"
#include "materials/equation_of_state.h"
#include "materials/material_manager.h"
#include "solvers/convective_term_contributions/convective_term_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/riemann_solver_setup.h"
#include "solvers/state_reconstruction/state_reconstruction_setup.h"
#include "user_specifications/compile_time_constants.h"
#include "utilities/helper_functions.h"

/**
 * @brief Discretization of the convective term solver using a finite-volume
 * procedure.
 */
class FiniteVolumeScheme : public ConvectiveTermSolver<FiniteVolumeScheme> {

  friend ConvectiveTermSolver;

  using RiemannSolverConcretization = RiemannSolverSetup::Concretize<
      FiniteVolumeSettings::riemann_solver>::type;
  using StateReconstructionConcretization =
      StateReconstructionSetup::Concretize<state_reconstruction_type>::type;
  RiemannSolverConcretization const riemann_solver_;
  StateReconstructionConcretization const state_reconstruction_;

  template <Direction DIR>
  void ComputeFluxes(
      std::pair<MaterialName const, Block> const &mat_block,
      double (&fluxes)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (&u_hllc)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double const (
          &Roe_eigenvectors_left)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                                 [MF::ANOE()][MF::ANOE()],
      double const (
          &Roe_eigenvectors_right)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                                  [MF::ANOE()][MF::ANOE()],
      double const cell_size) const;

  void UpdateImplementation(
      std::pair<MaterialName const, Block> const &mat_block,
      double const cell_size,
      double (
          &fluxes_x)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (
          &fluxes_y)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (
          &fluxes_z)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (
          &volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()]) const;

public:
  FiniteVolumeScheme() = delete;
  explicit FiniteVolumeScheme(
      MaterialManager const &material_manager,
      EigenDecomposition const &eigendecomposition_calculator);
  ~FiniteVolumeScheme() = default;
  FiniteVolumeScheme(FiniteVolumeScheme const &) = delete;
  FiniteVolumeScheme &operator=(FiniteVolumeScheme const &) = delete;
  FiniteVolumeScheme(FiniteVolumeScheme &&) = delete;
  FiniteVolumeScheme &operator=(FiniteVolumeScheme &&) = delete;
};

#endif // FINITE_VOLUME_SCHEME_H
