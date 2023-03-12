//===------------------- interface_term_solver.h --------------------------===//
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
#ifndef INTERFACE_TERM_SOLVER_H
#define INTERFACE_TERM_SOLVER_H

#include "interface_interaction/interface_term_contributions/heat_exchange_fluxes.h"
#include "interface_interaction/interface_term_contributions/interface_stress_tensor_fluxes.h"
#include "materials/material_manager.h"
#include "stencils/spatial_derivative_stencils/central_difference.h"
#include "stencils/spatial_derivative_stencils/fourth_order_central_difference.h"
#include "topology/node.h"
#include "user_specifications/numerical_setup.h"

#include "levelset/geometry/geometry_calculator_setup.h"

using GeometryCalculatorConcretization =
    GeometryCalculatorSetup::Concretize<geometry_calculator>::type;

/**
 * The InterfaceTermSolver class handles the solution of interface pressure,
 * velocity, and exchange terms.
 */
class InterfaceTermSolver {

private:
  static constexpr double one_twelfth_ = 1.0 / 12.0;

  MaterialManager const &material_manager_;
  GeometryCalculatorConcretization const geometry_calculator_;
  InterfaceStressTensorFluxes const interface_stress_tensor_fluxes_;
  HeatExchangeFluxes const heat_exchange_fluxes_;

  void FillInterfaceNormalVelocityBuffer(
      Node const &node,
      double (
          &u_interface_normal_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3]) const;

  void FillDeltaApertureBuffer(
      Node const &node,
      double (&delta_aperture_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3]) const;

public:
  InterfaceTermSolver() = delete;
  explicit InterfaceTermSolver(MaterialManager const &material_manager);
  ~InterfaceTermSolver() = default;
  InterfaceTermSolver(InterfaceTermSolver const &) = delete;
  InterfaceTermSolver &operator=(InterfaceTermSolver const &) = delete;
  InterfaceTermSolver(InterfaceTermSolver &&) = delete;
  InterfaceTermSolver &operator=(InterfaceTermSolver &&) = delete;

  void SolveInterfaceInteraction(Node &node) const;
  void
  WeightFaceFluxes(Node const &node, MaterialName const material,
                   double (&face_fluxes_x)[MF::ANOE()][CC::ICX() + 1]
                                          [CC::ICY() + 1][CC::ICZ() + 1],
                   double (&face_fluxes_y)[MF::ANOE()][CC::ICX() + 1]
                                          [CC::ICY() + 1][CC::ICZ() + 1],
                   double (&face_fluxes_z)[MF::ANOE()][CC::ICX() + 1]
                                          [CC::ICY() + 1][CC::ICZ() + 1]) const;

  void WeightVolumeForces(Node const &node, MaterialName const material,
                          double (&volume_forces)[MF::ANOE()][CC::ICX()]
                                                 [CC::ICY()][CC::ICZ()]) const;
};

#endif // INTERFACE_TERM_SOLVER_H
