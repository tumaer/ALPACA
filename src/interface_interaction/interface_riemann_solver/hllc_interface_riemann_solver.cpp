//===---------------- hllc_interface_riemann_solver.cpp -------------------===//
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
#include "hllc_interface_riemann_solver.h"
#include <cmath>
#include <limits>

/**
 * @brief Constructor for the HllcInterfaceRiemannSolver.
 * @param material_manager See base class.
 */
HllcInterfaceRiemannSolver::HllcInterfaceRiemannSolver(
    MaterialManager const &material_manager)
    : InterfaceRiemannSolver(material_manager) {
  // Empty besides call of base class constructor.
}

/**
 * @brief Solves the Riemann problem at the interface with a generalized HLLC
 * approach. See \cite Hu2008 for details.
 * @note DO NOT USE with WATERLIKE EOS.
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
HllcInterfaceRiemannSolver::SolveInterfaceRiemannProblemImplementation(
    double const rho_left, double const p_left,
    double const velocity_normal_left, MaterialName const material_left,
    double const rho_right, double const p_right,
    double const velocity_normal_right, MaterialName const material_right,
    double const) const {
  // obtain parameters of the left (negative) material
  double const one_rho_left =
      1.0 / std::max(rho_left, std::numeric_limits<double>::epsilon());
  double const c_left = material_manager_.GetMaterial(material_left)
                            .GetEquationOfState()
                            .SpeedOfSound(rho_left, p_left);
  double const gruneisen_left = material_manager_.GetMaterial(material_left)
                                    .GetEquationOfState()
                                    .Gruneisen(rho_left);
  double const psi_left = material_manager_.GetMaterial(material_left)
                              .GetEquationOfState()
                              .Psi(p_left, one_rho_left);
  double const gamma_left =
      material_manager_.GetMaterial(material_left).GetEquationOfState().Gamma();
  double const b_left =
      material_manager_.GetMaterial(material_left).GetEquationOfState().B();

  // obtain parameters of the right (positive) material
  double const one_rho_right =
      1.0 / std::max(rho_right, std::numeric_limits<double>::epsilon());
  double const c_right = material_manager_.GetMaterial(material_right)
                             .GetEquationOfState()
                             .SpeedOfSound(rho_right, p_right);
  double const gruneisen_right = material_manager_.GetMaterial(material_right)
                                     .GetEquationOfState()
                                     .Gruneisen(rho_right);
  double const psi_right = material_manager_.GetMaterial(material_right)
                               .GetEquationOfState()
                               .Psi(p_right, one_rho_right);
  double const gamma_right = material_manager_.GetMaterial(material_right)
                                 .GetEquationOfState()
                                 .Gamma();
  double const b_right =
      material_manager_.GetMaterial(material_right).GetEquationOfState().B();

  // compute expensive and frequently used temporaries
  double const sqrt_rho_left =
      std::sqrt(std::max(rho_left, std::numeric_limits<double>::epsilon()));
  double const sqrt_rho_right =
      std::sqrt(std::max(rho_right, std::numeric_limits<double>::epsilon()));
  double const rho_div = 1.0 / (sqrt_rho_left + sqrt_rho_right);

  // compute Roe averages
  double const density_roe_ave = sqrt_rho_left * sqrt_rho_right;
  double const velocity_roe_ave = ((velocity_normal_left * sqrt_rho_left) +
                                   (velocity_normal_right * sqrt_rho_right)) *
                                  rho_div;
  double const psi_roe_ave =
      ((psi_left * sqrt_rho_left) + (psi_right * sqrt_rho_right)) * rho_div;
  double const gruneisen_roe_ave =
      ((gruneisen_left * sqrt_rho_left) + (gruneisen_right * sqrt_rho_right)) *
      rho_div;

  // compute Roe averaged speed of sound
  double const tmp = (velocity_normal_right - velocity_normal_left) *
                     (velocity_normal_right - velocity_normal_left);
  double const pressure_over_density_roe_ave =
      ((p_left * one_rho_left) * sqrt_rho_left +
       (p_right * one_rho_right) * sqrt_rho_right) *
          rho_div +
      0.5 * density_roe_ave * (rho_div * rho_div) * tmp;
  double const c_roe_ave = std::sqrt(
      psi_roe_ave + gruneisen_roe_ave * pressure_over_density_roe_ave);

  // compute signal speeds
  double const signal_speed_left =
      std::min(velocity_normal_left - c_left, velocity_roe_ave - c_roe_ave);
  double const signal_speed_right =
      std::max(velocity_normal_right + c_right, velocity_roe_ave + c_roe_ave);

  // compute interface velocity
  double const interface_velocity =
      (rho_right * velocity_normal_right *
           (signal_speed_right - velocity_normal_right) +
       rho_left * velocity_normal_left *
           (velocity_normal_left - signal_speed_left) +
       (p_left - p_right)) /
      (rho_right * (signal_speed_right - velocity_normal_right) +
       rho_left * (velocity_normal_left - signal_speed_left));

  // compute interface pressure
  double const alpha = (interface_velocity - signal_speed_left) /
                       (signal_speed_right - signal_speed_left);
  double const beta = (signal_speed_right - interface_velocity) /
                      (signal_speed_right - signal_speed_left);
  double const rho_left_star =
      rho_left * (signal_speed_left - velocity_normal_left) /
      (signal_speed_left - interface_velocity); // incorrect in paper
  double const rho_right_star =
      rho_right * (signal_speed_right - velocity_normal_right) /
      (signal_speed_right - interface_velocity); // incorrect in paper
  double const energy_left =
      material_manager_.GetMaterial(material_left)
          .GetEquationOfState()
          .Energy(rho_left, velocity_normal_left, 0.0, 0.0, p_left);
  double const energy_right =
      material_manager_.GetMaterial(material_right)
          .GetEquationOfState()
          .Energy(rho_right, velocity_normal_right, 0.0, 0.0, p_right);
  double const energy_star =
      1.0 / (signal_speed_right - signal_speed_left) *
          ((energy_right * signal_speed_right -
            energy_left * signal_speed_left) +
           ((energy_left + p_left) * velocity_normal_left -
            (energy_right + p_right) * velocity_normal_right)) -
      0.5 * (alpha * rho_left_star + beta * rho_right_star) *
          (interface_velocity * interface_velocity);
  double const f_l = -gamma_left * b_left;   // only valid for stiffened EOS
  double const f_r = -gamma_right * b_right; // only valid for stiffened EOS
  double const interface_pressure_positive =
      (gruneisen_left * gruneisen_right) /
      (beta * gruneisen_left + alpha * gruneisen_right) *
      (energy_star +
       (alpha * (f_l / gruneisen_left) + beta * (f_r / gruneisen_right)));

  double const interface_pressure_negative = interface_pressure_positive;

  return {interface_velocity, interface_pressure_positive,
          interface_pressure_negative};
}
