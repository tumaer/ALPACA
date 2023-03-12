//===------------ linearized_interface_riemann_solver.cpp -----------------===//
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
#include "linearized_interface_riemann_solver.h"
#include <limits>

/**
 * @brief Constructor for the LinearizedInterfaceRiemannSolver.
 * @param material_manager See base class.
 */
LinearizedInterfaceRiemannSolver::LinearizedInterfaceRiemannSolver(
    MaterialManager const &material_manager)
    : InterfaceRiemannSolver(material_manager) {
  // Empty besides call of base class constructor.
}

/**
 * @brief Solves the Riemann problem at the interface linearized. See \cite
 * Saurel2003 for details.
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
 * interface_velocity, interface_pressure_positive, interface_pressure_negative.
 */
std::array<double, 3>
LinearizedInterfaceRiemannSolver::SolveInterfaceRiemannProblemImplementation(
    double const rho_left, double const p_left,
    double const velocity_normal_left, MaterialName const material_left,
    double const rho_right, double const p_right,
    double const velocity_normal_right, MaterialName const material_right,
    double const delta_p) const {
  double const c_left = material_manager_.GetMaterial(material_left)
                            .GetEquationOfState()
                            .SpeedOfSound(rho_left, p_left);
  double const impedance_left = rho_left * c_left;

  double const c_right = material_manager_.GetMaterial(material_right)
                             .GetEquationOfState()
                             .SpeedOfSound(rho_right, p_right);
  double const impedance_right = rho_right * c_right;

  double const inverse_impedance_sum =
      1.0 / std::max((impedance_left + impedance_right),
                     std::numeric_limits<double>::epsilon());

  double const interface_velocity =
      (impedance_left * velocity_normal_left +
       impedance_right * velocity_normal_right + p_left - p_right - delta_p) *
      inverse_impedance_sum;

  double const interface_pressure_positive =
      (impedance_left * p_right + impedance_right * (p_left - delta_p) +
       impedance_left * impedance_right *
           (velocity_normal_left - velocity_normal_right)) *
      inverse_impedance_sum;

  double const interface_pressure_negative =
      (impedance_left * (p_right + delta_p) + impedance_right * p_left +
       impedance_left * impedance_right *
           (velocity_normal_left - velocity_normal_right)) *
      inverse_impedance_sum;

  return {interface_velocity, interface_pressure_positive,
          interface_pressure_negative};
}
