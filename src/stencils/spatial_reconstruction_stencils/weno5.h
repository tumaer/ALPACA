//===----------------------------- weno5.h --------------------------------===//
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
#ifndef WENO5_H
#define WENO5_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute
 * fluxes according to \cite Jiang1996.
 */
class WENO5 : public Stencil<WENO5> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

  // Coefficients for WENO5 scheme
  static constexpr double coef_smoothness_1_ = 13.0 / 12.0;
  static constexpr double coef_smoothness_2_ = 0.25;

  static constexpr double coef_smoothness_11_ = 1.0;
  static constexpr double coef_smoothness_12_ = -2.0;
  static constexpr double coef_smoothness_13_ = 1.0;
  static constexpr double coef_smoothness_14_ = 1.0;
  static constexpr double coef_smoothness_15_ = -4.0;
  static constexpr double coef_smoothness_16_ = 3.0;

  static constexpr double coef_smoothness_21_ = 1.0;
  static constexpr double coef_smoothness_22_ = -2.0;
  static constexpr double coef_smoothness_23_ = 1.0;
  static constexpr double coef_smoothness_24_ = 1.0;
  static constexpr double coef_smoothness_25_ = -1.0;

  static constexpr double coef_smoothness_31_ = 1.0;
  static constexpr double coef_smoothness_32_ = -2.0;
  static constexpr double coef_smoothness_33_ = 1.0;
  static constexpr double coef_smoothness_34_ = 3.0;
  static constexpr double coef_smoothness_35_ = -4.0;
  static constexpr double coef_smoothness_36_ = 1.0;

  static constexpr double coef_weights_1_ = 0.1;
  static constexpr double coef_weights_2_ = 0.6;
  static constexpr double coef_weights_3_ = 0.3;

  static constexpr double coef_stencils_1_ = 2.0 / 6.0;
  static constexpr double coef_stencils_2_ = -7.0 / 6.0;
  static constexpr double coef_stencils_3_ = 11.0 / 6.0;
  static constexpr double coef_stencils_4_ = -1.0 / 6.0;
  static constexpr double coef_stencils_5_ = 5.0 / 6.0;
  static constexpr double coef_stencils_6_ = 2.0 / 6.0;
  static constexpr double coef_stencils_7_ = 2.0 / 6.0;
  static constexpr double coef_stencils_8_ = 5.0 / 6.0;
  static constexpr double coef_stencils_9_ = -1.0 / 6.0;

  // Small values to avoid division by 0, but also to adjust dissipation.
  static constexpr double epsilon_weno5_ = 1.0e-6;

  // Number of cells required for upwind and downwind stencils, as well as
  // number of cells downstream of the cell
  static constexpr unsigned int stencil_size_ = 6;
  static constexpr unsigned int downstream_stencil_size_ = 2;

  /**
   * @brief Evaluates the stencil according to a WENO-5 scheme. Also See base
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
              2 * evaluation_properties[1]];
    double const v2 =
        array[downstream_stencil_size_ + evaluation_properties[0] -
              1 * evaluation_properties[1]];
    double const v3 =
        array[downstream_stencil_size_ + evaluation_properties[0]];
    double const v4 =
        array[downstream_stencil_size_ + evaluation_properties[0] +
              1 * evaluation_properties[1]];
    double const v5 =
        array[downstream_stencil_size_ + evaluation_properties[0] +
              2 * evaluation_properties[1]];

    // Compute smoothness indicators s_i
    double const s11 = coef_smoothness_11_ * v1 + coef_smoothness_12_ * v2 +
                       coef_smoothness_13_ * v3;
    double const s12 = coef_smoothness_14_ * v1 + coef_smoothness_15_ * v2 +
                       coef_smoothness_16_ * v3;

    double const s1 =
        coef_smoothness_1_ * s11 * s11 + coef_smoothness_2_ * s12 * s12;

    double const s21 = coef_smoothness_21_ * v2 + coef_smoothness_22_ * v3 +
                       coef_smoothness_23_ * v4;
    double const s22 = coef_smoothness_24_ * v2 + coef_smoothness_25_ * v4;

    double const s2 =
        coef_smoothness_1_ * s21 * s21 + coef_smoothness_2_ * s22 * s22;

    double const s31 = coef_smoothness_31_ * v3 + coef_smoothness_32_ * v4 +
                       coef_smoothness_33_ * v5;
    double const s32 = coef_smoothness_34_ * v3 + coef_smoothness_35_ * v4 +
                       coef_smoothness_36_ * v5;

    double const s3 =
        coef_smoothness_1_ * s31 * s31 + coef_smoothness_2_ * s32 * s32;

    // Compute weights
    // NOTE: The epsilon value is used here explicitly to avoid compiler
    // optimizations when the epsilon is added directly to s_i.
    //       This could lead to undesired behavior in case the values s1, ...,
    //       s3 are of similar magnitude. Then, it cannot guaranteed anymore
    //       that a division by zero is avoided.
    double const a1 =
        coef_weights_1_ / ((s1 + epsilon_weno5_) * (s1 + epsilon_weno5_));
    double const a2 =
        coef_weights_2_ / ((s2 + epsilon_weno5_) * (s2 + epsilon_weno5_));
    double const a3 =
        coef_weights_3_ / ((s3 + epsilon_weno5_) * (s3 + epsilon_weno5_));

    double const one_a_sum = 1.0 / (a1 + a2 + a3);

    double const w1 = a1 * one_a_sum;
    double const w2 = a2 * one_a_sum;
    double const w3 = a3 * one_a_sum;

    // Return weighted average
    return w1 * (coef_stencils_1_ * v1 + coef_stencils_2_ * v2 +
                 coef_stencils_3_ * v3) +
           w2 * (coef_stencils_4_ * v2 + coef_stencils_5_ * v3 +
                 coef_stencils_6_ * v4) +
           w3 * (coef_stencils_7_ * v3 + coef_stencils_8_ * v4 +
                 coef_stencils_9_ * v5);
  }

public:
  explicit constexpr WENO5() = default;
  ~WENO5() = default;
  WENO5(WENO5 const &) = delete;
  WENO5 &operator=(WENO5 const &) = delete;
  WENO5(WENO5 &&) = delete;
  WENO5 &operator=(WENO5 &&) = delete;
};

#endif // STENCIL_WENO5_H
