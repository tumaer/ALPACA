//===------------------------ space_solver.cpp ----------------------------===//
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
#include "space_solver.h"

#include "utilities/mathematical_functions.h"
#include <cmath>

/**
 * @brief Standard constructor using an already existing MaterialManager and the
 * user-defined gravity.
 * @param material_manager The material manager object provides the correct
 * equation of state for the given material.
 * @param gravity Three-dimensional array holding the gravitational pull in x-,
 * y-, z-direction.
 */
SpaceSolver::SpaceSolver(MaterialManager const &material_manager,
                         std::array<double, 3> const gravity)
    : eigendecomposition_calculator_(material_manager),
      convective_term_solver_(material_manager, eigendecomposition_calculator_),
      source_term_solver_(material_manager, gravity),
      interface_term_solver_(material_manager),
      material_manager_(material_manager), levelset_advector_() {
  /* Empty besides initializer list*/
}

/**
 * @brief Computes right side of the underlying system of equations (including
 * source terms).
 * @param node The node under consideration.
 */
void SpaceSolver::UpdateFluxes(Node &node) const {

  double face_fluxes_x[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1];
  double face_fluxes_y[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1];
  double face_fluxes_z[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1];

  double volume_forces[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];

  for (auto &phase : node.GetPhases()) {
    for (Equation const eq : MF::ASOE()) {
      double(&rhs_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          phase.second.GetRightHandSideBuffer(eq);
      for (unsigned int i = 0; i < CC::TCX(); ++i) {
        for (unsigned int j = 0; j < CC::TCY(); ++j) {
          for (unsigned int k = 0; k < CC::TCZ(); ++k) {
            rhs_buffer[i][j][k] = 0.0;
          } // k
        }   // j
      }     // i
    }       // equation
  }         // phases

  // For multi-phase + Lmax nodes (which have a levelset)
  if (node.HasLevelset()) {
    // Solve interface Riemann problem to obtain interface velocity and
    // interface exchange terms
    interface_term_solver_.SolveInterfaceInteraction(node);
  }

  double const cell_size = node.GetCellSize();
  double const one_cell_size = 1.0 / cell_size;

  for (auto &phase : node.GetPhases()) {
    if constexpr (CC::SolidBoundaryActive()) {
      if (material_manager_.IsSolidBoundary(phase.first))
        continue;
    }

    // RHS-buffers have to be reset for each phase!
    for (unsigned int e = 0; e < MF::ANOE(); ++e) {
      for (unsigned int i = 0; i < CC::ICX() + 1; ++i) {
        for (unsigned int j = 0; j < CC::ICY() + 1; ++j) {
          for (unsigned int k = 0; k < CC::ICZ() + 1; ++k) {
            face_fluxes_x[e][i][j][k] = 0.0;
            face_fluxes_y[e][i][j][k] = 0.0;
            face_fluxes_z[e][i][j][k] = 0.0;
          } // k
        }   // j
      }     // i
    }       // equation

    for (unsigned int e = 0; e < MF::ANOE(); ++e) {
      for (unsigned int i = 0; i < CC::ICX(); ++i) {
        for (unsigned int j = 0; j < CC::ICY(); ++j) {
          for (unsigned int k = 0; k < CC::ICZ(); ++k) {
            volume_forces[e][i][j][k] = 0.0;
          } // k
        }   // j
      }     // i
    }       // equation

    // Determine cell face fluxes unsing a Riemann solver
    if constexpr (CC::InviscidExchangeActive()) {
      convective_term_solver_.UpdateConvectiveFluxes(
          phase, cell_size, face_fluxes_x, face_fluxes_y, face_fluxes_z,
          volume_forces);
    }

    // Determine source terms
    source_term_solver_.Sources(
        phase, cell_size, std::get<0>(node.GetBlockCoordinates()),
        face_fluxes_x, face_fluxes_y, face_fluxes_z, volume_forces);

    if (node.HasLevelset()) {
      interface_term_solver_.WeightFaceFluxes(node, phase.first, face_fluxes_x,
                                              face_fluxes_y, face_fluxes_z);
      interface_term_solver_.WeightVolumeForces(node, phase.first,
                                                volume_forces);
    }

    // update cells due to fluxes
    for (Equation const eq : MF::ASOE()) {
      auto const e = ETI(eq);
      double(&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          phase.second.GetRightHandSideBuffer(eq);
      for (unsigned int i = 0; i < CC::ICX(); ++i) {
        for (unsigned int j = 0; j < CC::ICY(); ++j) {
          for (unsigned int k = 0; k < CC::ICZ(); ++k) {
            // We need to add here because the RHS buffer might already be
            // filled in cut cells for single-level-set simulations
            cells[i + CC::FICX()][j + CC::FICY()][k + CC::FICZ()] +=
                volume_forces[e][i][j][k] +
                DimensionAwareConsistencyManagedSum(
                    face_fluxes_x[e][i][j + 1][k + 1] -
                        face_fluxes_x[e][i + 1][j + 1][k + 1],
                    face_fluxes_y[e][i + 1][j][k + 1] -
                        face_fluxes_y[e][i + 1][j + 1][k + 1],
                    face_fluxes_z[e][i + 1][j + 1][k] -
                        face_fluxes_z[e][i + 1][j + 1][k + 1]) *
                    one_cell_size;
          } // k
        }   // j
      }     // i
    }       // equation

    // save boundary fluxes for correction at jump boundary conditions
    double(&boundary_fluxes_west)[MF::ANOE()][CC::ICY()][CC::ICZ()] =
        phase.second.GetBoundaryJumpFluxes(BoundaryLocation::West);
    double(&boundary_fluxes_east)[MF::ANOE()][CC::ICY()][CC::ICZ()] =
        phase.second.GetBoundaryJumpFluxes(BoundaryLocation::East);
    double(&boundary_fluxes_south)[MF::ANOE()][CC::ICY()][CC::ICZ()] =
        phase.second.GetBoundaryJumpFluxes(BoundaryLocation::South);
    double(&boundary_fluxes_north)[MF::ANOE()][CC::ICY()][CC::ICZ()] =
        phase.second.GetBoundaryJumpFluxes(BoundaryLocation::North);
    double(&boundary_fluxes_bottom)[MF::ANOE()][CC::ICY()][CC::ICZ()] =
        phase.second.GetBoundaryJumpFluxes(BoundaryLocation::Bottom);
    double(&boundary_fluxes_top)[MF::ANOE()][CC::ICY()][CC::ICZ()] =
        phase.second.GetBoundaryJumpFluxes(BoundaryLocation::Top);
    for (unsigned int e = 0; e < MF::ANOE(); ++e) {
      for (unsigned int i = 0; i < CC::ICY(); ++i) {
        for (unsigned int j = 0; j < CC::ICZ(); ++j) {
          boundary_fluxes_west[e][i][j] += face_fluxes_x[e][0][i + 1][j + 1];
          boundary_fluxes_east[e][i][j] +=
              face_fluxes_x[e][CC::ICX()][i + 1][j + 1];
          boundary_fluxes_south[e][i][j] += face_fluxes_y[e][i + 1][0][j + 1];
          boundary_fluxes_north[e][i][j] +=
              face_fluxes_y[e][i + 1][CC::ICY()][j + 1];
          boundary_fluxes_bottom[e][i][j] += face_fluxes_z[e][i + 1][j + 1][0];
          boundary_fluxes_top[e][i][j] +=
              face_fluxes_z[e][i + 1][j + 1][CC::ICZ()];
        }
      }
    }
  } // phases
}

/**
 * @brief Computes right-hand side of the level set advection equation.
 * @param node The node under consideration.
 */
void SpaceSolver::UpdateLevelsetFluxes(Node &node) const {

  double(&levelset_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock().GetRightHandSideBuffer(
          InterfaceDescription::Levelset);
  for (unsigned int i = 0; i < CC::TCX(); ++i) {
    for (unsigned int j = 0; j < CC::TCY(); ++j) {
      for (unsigned int k = 0; k < CC::TCZ(); ++k) {
        levelset_rhs[i][j][k] = 0.0;
      } // k
    }   // j
  }     // i

  // compute levelset rhs (from levelset and interface velocity)
  levelset_advector_.Advect(node);
}

/**
 * @brief Computes the maximum eigenvalue in the given node.
 * @param block The block for which the eigenvalues are to be computed.
 * @param eigenvalues Indirect return parameter.
 */
void SpaceSolver::ComputeMaxEigenvaluesForPhase(
    std::pair<MaterialName const, Block> const &mat_block,
    double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()]) const {
  eigendecomposition_calculator_.ComputeMaxEigenvaluesOnBlock(mat_block,
                                                              eigenvalues);
}

/**
 * @brief Stores the given (GLF) eigenvalues for later usage in a globally valid
 * location.
 * @param eigenvalues Values to be stored.
 */
void SpaceSolver::SetFluxFunctionGlobalEigenvalues(
    double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()]) const {
  eigendecomposition_calculator_.SetGlobalEigenvalues(eigenvalues);
}
