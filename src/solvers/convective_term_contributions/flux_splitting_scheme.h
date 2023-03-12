//===--------------------- flux_splitting_scheme.h ------------------------===//
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
#ifndef FLUX_SPLITTING_SCHEME_H
#define FLUX_SPLITTING_SCHEME_H

#include "block_definitions/block.h"
#include "enums/direction_definition.h"
#include "materials/equation_of_state.h"
#include "materials/material_manager.h"
#include "solvers/convective_term_contributions/convective_term_solver.h"
#include "user_specifications/compile_time_constants.h"
#include "utilities/helper_functions.h"

/**
 * @brief Discretization of the convective term solver using a flux-splitting
 * procedure.
 */
class FluxSplittingScheme : public ConvectiveTermSolver<FluxSplittingScheme> {

  friend ConvectiveTermSolver;

  template <Direction DIR>
  void ComputeFluxes(
      Block const &b,
      double (&fluxes)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (&advection)[MF::ANOE()][CC::TCX()][CC::TCY()][CC::TCZ()],
      double const cell_size,
      double (&roe_eigenvectors_left)[CC::ICX() + 1][CC::ICY() + 1]
                                     [CC::ICZ() + 1][MF::ANOE()][MF::ANOE()],
      double (&roe_eigenvectors_right)[CC::ICX() + 1][CC::ICY() + 1]
                                      [CC::ICZ() + 1][MF::ANOE()][MF::ANOE()],
      double (&fluxfunction_wavespeed)[CC::ICX() + 1][CC::ICY() + 1]
                                      [CC::ICZ() + 1][MF::ANOE()]) const;

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
  FluxSplittingScheme() = delete;
  explicit FluxSplittingScheme(
      MaterialManager const &material_manager,
      EigenDecomposition const &eigendecomposition_calculator);
  ~FluxSplittingScheme() = default;
  FluxSplittingScheme(FluxSplittingScheme const &) = delete;
  FluxSplittingScheme &operator=(FluxSplittingScheme const &) = delete;
  FluxSplittingScheme(FluxSplittingScheme &&) = delete;
  FluxSplittingScheme &operator=(FluxSplittingScheme &&) = delete;
};

#endif // FLUX_SPLITTING_SCHEME_H
