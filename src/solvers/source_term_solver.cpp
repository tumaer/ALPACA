//===--------------------- source_term_solver.cpp -------------------------===//
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
#include "source_term_solver.h"

/**
 * @brief Standard constructor using an already existing MaterialManager and the
 * user-defined gravity.
 * @param material_manager The material manager object provides the correct
 * equation of state for the given material.
 * @param gravity Three-dimensional array holding the gravitational pull in x-,
 * y-, z-direction.
 */
SourceTermSolver::SourceTermSolver(MaterialManager const &material_manager,
                                   std::array<double, 3> const gravity)
    : gravity_(gravity), viscous_fluxes_(material_manager),
      heat_fluxes_(material_manager), axisymmetric_fluxes_(),
      axisymmetric_viscous_volume_forces_(material_manager) {
  /* Empty besides initializer list*/
}

/**
 * @brief Computes additions to the right-hand side solution due to the present
 * source terms.
 * @param mat_block The phase with its material identifier.
 * @param cell_size The cell size of the node.
 * @param node_origin_x The coordinate of the node origin in x-direction.
 * @param face_fluxes_x, face_fluxes_y, face_fluxes_z Fluxes across the cell
 * face.
 * @param volume_forces The volume forces acting at the cell center of the node.
 */
void SourceTermSolver::Sources(
    std::pair<MaterialName const, Block> const &mat_block,
    double const cell_size, double const node_origin_x,
    double (&face_fluxes_x)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                           [CC::ICZ() + 1],
    double (&face_fluxes_y)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                           [CC::ICZ() + 1],
    double (&face_fluxes_z)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                           [CC::ICZ() + 1],
    double (
        &volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()]) const {

  // compute dissipative fluxes
  if constexpr (CC::ViscosityIsActive()) {
    viscous_fluxes_.ComputeFluxes(mat_block, face_fluxes_x, face_fluxes_y,
                                  face_fluxes_z, cell_size);
  }

  // compute changes due to gravity
  if constexpr (CC::GravityIsActive()) {
    gravity_.ComputeForces(mat_block.second, volume_forces);
  }

  // compute terms for axisymmetric simulations
  if constexpr (CC::Axisymmetric()) {
    axisymmetric_fluxes_.ComputeAxisymmetricContributions(
        mat_block.second, volume_forces, cell_size, node_origin_x);
  }

  if constexpr (CC::ViscosityIsActive() && CC::Axisymmetric()) {
    axisymmetric_viscous_volume_forces_.ComputeForces(mat_block, volume_forces,
                                                      cell_size, node_origin_x);
  }

  // Compute terms for heat exchange
  if constexpr (CC::HeatConductionActive() &&
                MF::IsEquationActive(Equation::Energy) &&
                MF::IsPrimeStateActive(PrimeState::Temperature)) {
    heat_fluxes_.ComputeFluxes(mat_block, face_fluxes_x, face_fluxes_y,
                               face_fluxes_z, cell_size);
  }
}
