//===-------------------- interface_riemann_solver.h ----------------------===//
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
#ifndef INTERFACE_RIEMANN_SOLVER_H
#define INTERFACE_RIEMANN_SOLVER_H

#include "materials/material_manager.h"

/**
 * @brief The class InterfaceRiemannSolver provides functionality to solve the
 * Riemann problem at the interface.
 * @tparam DerivedInterfaceRiemannSolver The template for the derived classes.
 * This is necessary for the CRTP.
 */
template <typename DerivedInterfaceRiemannSolver> class InterfaceRiemannSolver {

  friend DerivedInterfaceRiemannSolver;

protected:
  MaterialManager const &material_manager_;

  /**
   * @brief Constructor for the InterfaceRiemannSolver class.
   * @param material_manager Contains information about the materials present in
   * the simulation.
   */
  explicit InterfaceRiemannSolver(MaterialManager const &material_manager)
      : material_manager_(material_manager) {
    // Empty besides initializer list.
  }

public:
  InterfaceRiemannSolver() = delete;
  ~InterfaceRiemannSolver() = default;
  InterfaceRiemannSolver(InterfaceRiemannSolver const &) = delete;
  InterfaceRiemannSolver &operator=(InterfaceRiemannSolver const &) = delete;
  InterfaceRiemannSolver(InterfaceRiemannSolver &&) = delete;
  InterfaceRiemannSolver &operator=(InterfaceRiemannSolver &&) = delete;

  /**
   * @brief Solves the Riemann problem at the interface linearized.
   * @param rho_left Density of the left fluid.
   * @param p_left Pressure of the left fluid.
   * @param velocity_normal_left Velocity normal to the interface of the left
   * fluid.
   * @param material_left Material of the left fluid.
   * @param rho_right Density of the right fluid.
   * @param p_right Pressure of the right fluid.
   * @param velocity_normal_right Velocity normal to the interface of the right
   * fluid.
   * @param material_right Material of the right fluid.
   * @param delta_p Pressure jump due to capillarity.
   * @return An array that contains following information in the given order:
   * interface_velocity, interface_pressure_positive,
   * interface_pressure_negative.
   */
  std::array<double, 3> SolveInterfaceRiemannProblem(
      double const &rho_left, double const &p_left,
      double const &velocity_normal_left, MaterialName const &material_left,
      double const &rho_right, double const &p_right,
      double const &velocity_normal_right, MaterialName const &material_right,
      double const &delta_p) const {
    if constexpr (CC::SolidBoundaryActive()) {
      if (material_manager_.IsSolidBoundary(material_left)) {
        return {velocity_normal_left, p_right, 0.0};
      }
      if (material_manager_.IsSolidBoundary(material_right)) {
        return {velocity_normal_right, 0.0, p_left};
      }
    }
    return static_cast<DerivedInterfaceRiemannSolver const &>(*this)
        .SolveInterfaceRiemannProblemImplementation(
            rho_left, p_left, velocity_normal_left, material_left, rho_right,
            p_right, velocity_normal_right, material_right, delta_p);
  }
};

#endif // INTERFACE_RIEMANN_SOLVER_H
