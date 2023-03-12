//===----------------------- source_term_solver.h -------------------------===//
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
#ifndef SOURCE_TERM_SOLVER_H
#define SOURCE_TERM_SOLVER_H

#include <array>

#include "block_definitions/block.h"
#include "materials/material_manager.h"
#include "solvers/source_term_contributions/axisymmetric_fluxes.h"
#include "solvers/source_term_contributions/axisymmetric_viscous_volume_forces.h"
#include "solvers/source_term_contributions/gravitational_force.h"
#include "solvers/source_term_contributions/heat_fluxes.h"
#include "solvers/source_term_contributions/viscous_fluxes.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief The SourceTermSolver class adds contributions due to source terms. The
 * solution of the right side of the Euler equations.
 */
class SourceTermSolver {

  GravitationalForce const gravity_;
  ViscousFluxes const viscous_fluxes_;
  HeatFluxes const heat_fluxes_;
  AxisymmetricFluxes const axisymmetric_fluxes_;
  AxisymmetricViscousVolumeForces const axisymmetric_viscous_volume_forces_;

public:
  SourceTermSolver() = delete;
  explicit SourceTermSolver(MaterialManager const &material_manager,
                            std::array<double, 3> const gravity);
  ~SourceTermSolver() = default;
  SourceTermSolver(SourceTermSolver const &) = delete;
  SourceTermSolver &operator=(SourceTermSolver const &) = delete;
  SourceTermSolver(SourceTermSolver &&) = delete;
  SourceTermSolver &operator=(SourceTermSolver &&) = delete;

  void Sources(std::pair<MaterialName const, Block> const &mat_block,
               double const cell_size, double const node_origin_x,
               double (&face_fluxes_x)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                                      [CC::ICZ() + 1],
               double (&face_fluxes_y)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                                      [CC::ICZ() + 1],
               double (&face_fluxes_z)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                                      [CC::ICZ() + 1],
               double (&volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()]
                                      [CC::ICZ()]) const;
};

#endif // SOURCE_TERM_SOLVER_H
