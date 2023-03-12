//===----------------------- runge_kutta_2_TVD.h --------------------------===//
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
#ifndef RUNGE_KUTTA_2_TVD_H
#define RUNGE_KUTTA_2_TVD_H

#include "time_integrator.h"

/**
 * @brief The RungeKutta2TVD class integrates in time using a total variation
 * diminishing two-step Runge-Kutta method. On paper the equations are I)    u^*
 * =       u^n +           +       dt * f(u^n) II)   u^(n+1) = 0.5 * u^n + 0.5 *
 * u^* + 0.5 * dt * f(u^*) To reduce the memory footprint the equations are
 * reordered to allow the integration using two buffers only.
 */
class RungeKutta2TVD : public TimeIntegrator<RungeKutta2TVD> {

  friend TimeIntegrator;

  static constexpr unsigned int number_of_stages_ = 2;

  static constexpr std::array<double, number_of_stages_>
      timestep_multiplier_jump_conservatives_ = {
          0.5, // first stage
          0.5  // second stage
      };

  static constexpr std::array<double, number_of_stages_>
      timestep_multiplier_conservatives_ = {
          1.0, // first stage
          0.5  // second stage
      };

  static constexpr std::array<std::array<double, 2>, number_of_stages_ - 1>
      buffer_multiplier_ = {{
          {0.5, 0.5} // second stage
      }};

public:
  RungeKutta2TVD() = delete;
  ~RungeKutta2TVD() = default;
  RungeKutta2TVD(RungeKutta2TVD const &) = delete;
  RungeKutta2TVD &operator=(RungeKutta2TVD const &) = delete;
  RungeKutta2TVD(RungeKutta2TVD &&) = delete;
  RungeKutta2TVD &operator=(RungeKutta2TVD &&) = delete;

  /**
   * @brief Constructor.
   * @param start_time Time when the simulation should start.
   */
  explicit RungeKutta2TVD(double const start_time = 0.0)
      : TimeIntegrator(start_time) {}
};

#endif // RUNGE_KUTTA_2_TVD_H
