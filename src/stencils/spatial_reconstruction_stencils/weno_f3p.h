//===--------------------------- weno_f3p.h -------------------------------===//
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
#ifndef WENOF3P_H
#define WENOF3P_H

#include "stencils/stencil.h"
#include <cmath>

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute
 * fluxes according to WENOF3+ \cite Gande2020.
 */
class WENOF3P : public Stencil<WENOF3P> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

  // Coefficients for WENO3 scheme
  static constexpr double coef_smoothness_11_ = -1.0;
  static constexpr double coef_smoothness_12_ = 1.0;
  static constexpr double coef_smoothness_21_ = -1.0;
  static constexpr double coef_smoothness_22_ = 1.0;

  static constexpr double coef_weights_1_ = 1.0 / 3.0;
  static constexpr double coef_weights_2_ = 2.0 / 3.0;

  static constexpr double coef_stencils_1_ = -1.0 / 2.0;
  static constexpr double coef_stencils_2_ = 3.0 / 2.0;
  static constexpr double coef_stencils_3_ = 1.0 / 2.0;
  static constexpr double coef_stencils_4_ = 1.0 / 2.0;

  static constexpr double coef_beta3_dot_ = 1.0 / 12.0;

  static constexpr double lambda_constant = 1.0 / 3.0;

  // Number of cells required for upwind and downwind stencils, as well as
  // number of cells downstream of the cell
  static constexpr unsigned int stencil_size_ = 4;
  static constexpr unsigned int downstream_stencil_size_ = 1;

  /**
   * @brief Evaluates the stencil according to a WENO-3 scheme. Also See base
   * class.
   * @note Hotpath function.
   */
  constexpr double
  ApplyImplementation(std::array<double, stencil_size_> const &array,
                      std::array<int const, 2> const evaluation_properties,
                      double const cell_size) const {
    // Assign values to v_i to make it easier to read
    double const v1 =
        array[downstream_stencil_size_ + evaluation_properties[0] -
              evaluation_properties[1]];
    double const v2 =
        array[downstream_stencil_size_ + evaluation_properties[0]];
    double const v3 =
        array[downstream_stencil_size_ + evaluation_properties[0] +
              evaluation_properties[1]];

    // Compute smoothness indicators s_i
    double const s11 = coef_smoothness_11_ * v1 + coef_smoothness_12_ * v2;
    double const beta0 = s11 * s11;

    double const s21 = coef_smoothness_21_ * v2 + coef_smoothness_22_ * v3;
    double const beta1 = s21 * s21;

    double const beta3_dot =
        coef_beta3_dot_ * (((v1 + v3) - 2.0 * v2) * ((v1 + v3) - 2.0 * v2)) +
        0.25 * ((v1 - v3) * (v1 - v3));

    // Compute weights
    double const tau = std::abs(0.5 * (beta0 + beta1) - beta3_dot);
    double const lambda = std::pow(cell_size, lambda_constant);

    double const tmp = lambda / (tau + epsilon_);
    double const a0 =
        coef_weights_1_ * (1.0 + (tau + epsilon_) / (beta0 + epsilon_) +
                           tmp * (beta0 + epsilon_));
    double const a1 =
        coef_weights_2_ * (1.0 + (tau + epsilon_) / (beta1 + epsilon_) +
                           tmp * (beta1 + epsilon_));

    double const one_a_sum = 1.0 / (a0 + a1);

    double const w1 = a0 * one_a_sum;
    double const w2 = a1 * one_a_sum;

    // Return weighted average
    return w1 * (coef_stencils_1_ * v1 + coef_stencils_2_ * v2) +
           w2 * (coef_stencils_3_ * v2 + coef_stencils_4_ * v3);
  }

public:
  explicit constexpr WENOF3P() = default;
  ~WENOF3P() = default;
  WENOF3P(WENOF3P const &) = delete;
  WENOF3P &operator=(WENOF3P const &) = delete;
  WENOF3P(WENOF3P &&) = delete;
  WENOF3P &operator=(WENOF3P &&) = delete;
};

#endif // STENCIL_WENOF3P_H
