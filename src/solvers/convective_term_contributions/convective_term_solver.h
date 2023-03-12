//===--------------------- convective_term_solver.h -----------------------===//
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
#ifndef CONVECTIVE_TERM_SOLVER_H
#define CONVECTIVE_TERM_SOLVER_H

#include "block_definitions/block.h"
#include "materials/material_manager.h"
#include "solvers/eigendecomposition.h"

/**
 * @brief The ConvectiveTermSolver calculates the convective term using a
 * Riemann solver of choice with an interchangable stencil.
 */
template <typename DerivedConvectiveTermSolver> class ConvectiveTermSolver {

  friend DerivedConvectiveTermSolver;

  EigenDecomposition const &eigendecomposition_calculator_;
  MaterialManager const &material_manager_;

  explicit ConvectiveTermSolver(
      MaterialManager const &material_manager,
      EigenDecomposition const &eigendecomposition_calculator)
      : eigendecomposition_calculator_(eigendecomposition_calculator),
        material_manager_(material_manager) {
    // Empty constructor besides initializer list.
  }

public:
  ConvectiveTermSolver() = delete;
  explicit ConvectiveTermSolver(MaterialManager const &material_manager,
                                std::array<double, 3> gravity);
  ~ConvectiveTermSolver() = default;
  ConvectiveTermSolver(ConvectiveTermSolver const &) = delete;
  ConvectiveTermSolver &operator=(ConvectiveTermSolver const &) = delete;
  ConvectiveTermSolver(ConvectiveTermSolver &&) = delete;
  ConvectiveTermSolver &operator=(ConvectiveTermSolver &&) = delete;

  /**
   * @brief Determines the contribution of the convective term of the underlying
   * system of equations within the provided node.
   * @param mat_block Container holding the relevant fluid data to compute the
   * update.
   * @param cell_size The size of the cells in the block.
   * @param fluxes_x, fluxes_y, fluxes_z The fluxes over the cell faces as
   * computed by this Riemann solver. Indirect return parameter.
   */
  void UpdateConvectiveFluxes(
      std::pair<MaterialName const, Block> const &mat_block,
      double const cell_size,
      double (
          &fluxes_x)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (
          &fluxes_y)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (
          &fluxes_z)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
      double (
          &volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()]) const {
    static_cast<DerivedConvectiveTermSolver const &>(*this)
        .UpdateImplementation(mat_block, cell_size, fluxes_x, fluxes_y,
                              fluxes_z, volume_forces);
  }
};

#endif // CONVECTIVE_TERM_SOLVER_H
