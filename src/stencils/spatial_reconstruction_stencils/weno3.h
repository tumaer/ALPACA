//===---------------------------- weno3.h ---------------------------------===//
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
#ifndef WENO3_H
#define WENO3_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute
 * fluxes according to \cite Shu1999.
 */
class WENO3 : public Stencil<WENO3> {

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

  // Small values to avoid division by 0, but also to adjust dissipation.
  // Optimized according to F. Schranner (same as WENO5)
  static constexpr double epsilon_1_ = 1.0e-6;
  static constexpr double epsilon_2_ = 1.0e-15;

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
                      double const) const {
    // Assign values to v_i to make it easier to read
    double const v1 =
        array[downstream_stencil_size_ + evaluation_properties[0] -
              1 * evaluation_properties[1]];
    double const v2 =
        array[downstream_stencil_size_ + evaluation_properties[0]];
    double const v3 =
        array[downstream_stencil_size_ + evaluation_properties[0] +
              1 * evaluation_properties[1]];

    // Compute smoothness indicators s_i
    double const s11 = coef_smoothness_11_ * v1 + coef_smoothness_12_ * v2;
    double const s1 = s11 * s11;

    double const s21 = coef_smoothness_21_ * v2 + coef_smoothness_22_ * v3;
    double const s2 = s21 * s21;

    // Compute weights
    // NOTE: The epsilon value is used here explicitly to avoid compiler
    // optimizations when the epsilon is added directly to s_i.
    //       This could lead to undesired behavior in case the values s1 and s2
    //       are of similar magnitude. Then, it cannot guaranteed anymore that a
    //       division by zero is avoided.
    double const a1 = coef_weights_1_ / ((s1 + epsilon_) * (s1 + epsilon_));
    double const a2 = coef_weights_2_ / ((s2 + epsilon_) * (s2 + epsilon_));

    double const one_a_sum = 1.0 / (a1 + a2);

    double const w1 = a1 * one_a_sum;
    double const w2 = a2 * one_a_sum;

    // Return weighted average
    return w1 * (coef_stencils_1_ * v1 + coef_stencils_2_ * v2) +
           w2 * (coef_stencils_3_ * v2 + coef_stencils_4_ * v3);
  }

public:
  explicit constexpr WENO3() = default;
  ~WENO3() = default;
  WENO3(WENO3 const &) = delete;
  WENO3 &operator=(WENO3 const &) = delete;
  WENO3(WENO3 &&) = delete;
  WENO3 &operator=(WENO3 &&) = delete;
};

#endif // STENCIL_WENO3_H
