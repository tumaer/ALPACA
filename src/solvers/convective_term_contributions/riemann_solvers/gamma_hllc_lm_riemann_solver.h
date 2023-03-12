//===----------------- gamma_hllc_lm_riemann_solver.h ---------------------===//
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
#ifndef GAMMA_HLLC_LM_RIEMANN_SOLVER_H
#define GAMMA_HLLC_LM_RIEMANN_SOLVER_H

#include "block_definitions/block.h"
#include "enums/direction_definition.h"
#include "materials/equation_of_state.h"
#include "materials/equations_of_state/gamma_model_stiffened_gas.h"
#include "materials/material_manager.h"
#include "solvers/convective_term_contributions/riemann_solvers/hll_signal_speed_calculator.h"
#include "solvers/convective_term_contributions/riemann_solvers/riemann_solver.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Discretization of the Riemann solver using the HLLC-LM procedure
 * according to \cite Fleischmann2020b with extension to gamma states according
 * to \Saurel2009.
 */
class GammaHllcLMRiemannSolver
    : public RiemannSolver<GammaHllcLMRiemannSolver> {

  friend RiemannSolver;

  static constexpr double Ma_limit = 0.1;

  /**
   * @brief Computes an approximate solution of a Riemann problem using the
   * HLLC-LM procedure.
   * @param material Material information of the phase under consideration.
   * @param conservatives_left The initial left conservative states.
   * @param conservatives_right The initial right conservative states.
   * @param prime_state_left The initial left primitive states.
   * @param prime_state_right The initial right primitive states.
   * @tparam DIR Indicates which spatial direction is to be computed.
   * @note Hotpath function.
   */
  template <Direction DIR>
  std::tuple<std::array<double, MF::ANOE()>, double>
  SolveGammaRiemannProblemImplementation(
      MaterialName const,
      std::array<double, MF::ANOE()> const &conservatives_left,
      std::array<double, MF::ANOE()> const &conservatives_right,
      std::array<double, MF::ANOP()> const &prime_state_left,
      std::array<double, MF::ANOP()> const &prime_state_right) const {
    std::array<double, MF::ANOE()> q_star_left;
    std::array<double, MF::ANOE()> q_star_right;
    std::array<double, MF::ANOE()> flux_left;
    std::array<double, MF::ANOE()> flux_right;
    std::array<double, MF::ANOE()> fluxes;

    constexpr unsigned int principal_momentum_index = ETI(MF::AME()[DTI(DIR)]);
    constexpr unsigned int principal_velocity_index = PTI(MF::AV()[DTI(DIR)]);

    // Compute pressure, velocity and speed of sound for both cells for
    // reconstructed values
    double const pressure_left = prime_state_left[PTI(PrimeState::Pressure)];
    double const pressure_right = prime_state_right[PTI(PrimeState::Pressure)];

    double const one_density_left =
        1.0 / conservatives_left[ETI(Equation::Mass)];
    double const one_density_right =
        1.0 / conservatives_right[ETI(Equation::Mass)];
    double const velocity_left = prime_state_left[principal_velocity_index];
    double const velocity_right = prime_state_right[principal_velocity_index];

    double const speed_of_sound_left =
        GammaModelStiffenedGas::CalculateSpeedOfSound(
            conservatives_left[ETI(Equation::Mass)], pressure_left,
            prime_state_left[PTI(PrimeState::gamma)],
            prime_state_left[PTI(PrimeState::pi)]);
    double const speed_of_sound_right =
        GammaModelStiffenedGas::CalculateSpeedOfSound(
            conservatives_right[ETI(Equation::Mass)], pressure_right,
            prime_state_right[PTI(PrimeState::gamma)],
            prime_state_right[PTI(PrimeState::pi)]);

    // Calculation of signal speeds
    auto const [wave_speed_left_simple, wave_speed_right_simple] =
        CalculateSignalSpeed(conservatives_left[ETI(Equation::Mass)],
                             conservatives_right[ETI(Equation::Mass)],
                             velocity_left, velocity_right, pressure_left,
                             pressure_right, speed_of_sound_left,
                             speed_of_sound_right, 0.0);

    double const wave_speed_contact =
        ((pressure_right - pressure_left) +
         (conservatives_left[principal_momentum_index] *
              (wave_speed_left_simple - velocity_left) -
          conservatives_right[principal_momentum_index] *
              (wave_speed_right_simple - velocity_right))) /
        (conservatives_left[ETI(Equation::Mass)] *
             (wave_speed_left_simple - velocity_left) -
         conservatives_right[ETI(Equation::Mass)] *
             (wave_speed_right_simple - velocity_right));
    double const wave_speed_left = std::min(wave_speed_left_simple, 0.0);
    double const wave_speed_right = std::max(wave_speed_right_simple, 0.0);

    double const chi_star_left = (wave_speed_left_simple - velocity_left) /
                                 (wave_speed_left_simple - wave_speed_contact);
    double const chi_star_right =
        (wave_speed_right_simple - velocity_right) /
        (wave_speed_right_simple - wave_speed_contact);

    // Compute intermediate states (Toro 10.71 10.72 10.73) and calculate F(q)
    q_star_left[ETI(Equation::Mass)] =
        conservatives_left[ETI(Equation::Mass)] * chi_star_left;
    q_star_left[principal_momentum_index] =
        conservatives_left[ETI(Equation::Mass)] * chi_star_left *
        wave_speed_contact;
    q_star_left[ETI(Equation::Energy)] =
        conservatives_left[ETI(Equation::Mass)] * chi_star_left *
        (conservatives_left[ETI(Equation::Energy)] * one_density_left +
         (wave_speed_contact - velocity_left) *
             (wave_speed_contact +
              pressure_left / (conservatives_left[ETI(Equation::Mass)] *
                               (wave_speed_left - velocity_left))));
    q_star_left[ETI(Equation::Gamma)] =
        conservatives_left[ETI(Equation::Gamma)] * chi_star_left;
    q_star_left[ETI(Equation::Pi)] =
        conservatives_left[ETI(Equation::Pi)] * chi_star_left;
    q_star_right[ETI(Equation::Mass)] =
        conservatives_right[ETI(Equation::Mass)] * chi_star_right;
    q_star_right[principal_momentum_index] =
        conservatives_right[ETI(Equation::Mass)] * chi_star_right *
        wave_speed_contact;
    q_star_right[ETI(Equation::Energy)] =
        conservatives_right[ETI(Equation::Mass)] * chi_star_right *
        (conservatives_right[ETI(Equation::Energy)] * one_density_right +
         (wave_speed_contact - velocity_right) *
             (wave_speed_contact +
              pressure_right / (conservatives_right[ETI(Equation::Mass)] *
                                (wave_speed_right - velocity_right))));
    q_star_right[ETI(Equation::Gamma)] =
        conservatives_right[ETI(Equation::Gamma)] * chi_star_right;
    q_star_right[ETI(Equation::Pi)] =
        conservatives_right[ETI(Equation::Pi)] * chi_star_right;

    flux_left[ETI(Equation::Mass)] =
        conservatives_left[principal_momentum_index];
    flux_left[principal_momentum_index] =
        ((conservatives_left[principal_momentum_index] *
          conservatives_left[principal_momentum_index]) *
         one_density_left) +
        pressure_left;
    flux_left[ETI(Equation::Energy)] =
        velocity_left *
        (conservatives_left[ETI(Equation::Energy)] + pressure_left);
    flux_left[ETI(Equation::Gamma)] =
        velocity_left * conservatives_left[ETI(Equation::Gamma)];
    flux_left[ETI(Equation::Pi)] =
        velocity_left * conservatives_left[ETI(Equation::Pi)];
    flux_right[ETI(Equation::Mass)] =
        conservatives_right[principal_momentum_index];
    flux_right[principal_momentum_index] =
        ((conservatives_right[principal_momentum_index] *
          conservatives_right[principal_momentum_index]) *
         one_density_right) +
        pressure_right;
    flux_right[ETI(Equation::Energy)] =
        velocity_right *
        (conservatives_right[ETI(Equation::Energy)] + pressure_right);
    flux_right[ETI(Equation::Gamma)] =
        velocity_right * conservatives_right[ETI(Equation::Gamma)];
    flux_right[ETI(Equation::Pi)] =
        velocity_right * conservatives_right[ETI(Equation::Pi)];

    // minor momenta
    for (unsigned int d = 0; d < DTI(CC::DIM()) - 1; ++d) {
      // get the index of this minor momentum
      unsigned int const minor_momentum_index =
          ETI(MF::AME()[DTI(GetMinorDirection<DIR>(d))]);

      q_star_left[minor_momentum_index] =
          chi_star_left * conservatives_left[minor_momentum_index];
      q_star_right[minor_momentum_index] =
          chi_star_right * conservatives_right[minor_momentum_index];

      flux_left[minor_momentum_index] =
          velocity_left * conservatives_left[minor_momentum_index];
      flux_right[minor_momentum_index] =
          velocity_right * conservatives_right[minor_momentum_index];
    }

    double const sig_left = 0.5 * (1.0 + Signum(wave_speed_contact));
    double const sig_right = 0.5 * (1.0 - Signum(wave_speed_contact));
    double const u_hllc = wave_speed_contact;

    for (unsigned int n = 0; n < MF::ANOE(); ++n) {
      double const flux_star_left =
          flux_left[n] +
          (wave_speed_left * (q_star_left[n] - conservatives_left[n]));
      double const flux_star_right =
          flux_right[n] +
          (wave_speed_right * (q_star_right[n] - conservatives_right[n]));
      fluxes[n] = sig_left * flux_star_left + sig_right * flux_star_right;
    }

    fluxes[ETI(Equation::Gamma)] =
        (sig_left * conservatives_left[ETI(Equation::Gamma)] *
         wave_speed_contact) +
        (sig_right * conservatives_right[ETI(Equation::Gamma)] *
         wave_speed_contact);
    fluxes[ETI(Equation::Pi)] =
        (sig_left * conservatives_left[ETI(Equation::Pi)] *
         wave_speed_contact) +
        (sig_right * conservatives_right[ETI(Equation::Pi)] *
         wave_speed_contact);

    double const Ma = std::min(
        Ma_limit, std::max(std::abs(velocity_left / speed_of_sound_left),
                           std::abs(velocity_right / speed_of_sound_right)));
    double const activation = std::sin((Ma / Ma_limit) * M_PI_2);

    for (unsigned int n = 0; n < MF::ANOE() - 2; ++n) {
      if (wave_speed_left >= 0.0) {
        fluxes[n] = flux_left[n];
      } else if (wave_speed_right <= 0.0) {
        fluxes[n] = flux_right[n];
      } else {
        fluxes[n] = 0.5 * ((flux_left[n] + flux_right[n]) -
                           (((activation * wave_speed_right) *
                                 (conservatives_right[n] - q_star_right[n]) +
                             (activation * wave_speed_left) *
                                 (conservatives_left[n] - q_star_left[n])) +
                            std::abs(wave_speed_contact) *
                                (q_star_right[n] - q_star_left[n])));
      }
    }
    return std::make_tuple(fluxes, u_hllc);
  }

public:
  GammaHllcLMRiemannSolver() = delete;
  explicit GammaHllcLMRiemannSolver(
      MaterialManager const &material_manager,
      EigenDecomposition const &eigendecomposition_calculator)
      : RiemannSolver(material_manager, eigendecomposition_calculator) {}
  ~GammaHllcLMRiemannSolver() = default;
  GammaHllcLMRiemannSolver(GammaHllcLMRiemannSolver const &) = delete;
  GammaHllcLMRiemannSolver &
  operator=(GammaHllcLMRiemannSolver const &) = delete;
  GammaHllcLMRiemannSolver(GammaHllcLMRiemannSolver &&) = delete;
  GammaHllcLMRiemannSolver &operator=(GammaHllcLMRiemannSolver &&) = delete;
};

#endif // GAMMA_HLLC_LM_RIEMANN_SOLVER_H
