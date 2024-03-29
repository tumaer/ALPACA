//===-------------------------- weno-ao53.h -------------------------------===//
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
#ifndef WENOAO53_H
#define WENOAO53_H

#include "stencils/stencil.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute
 * fluxes according to \cite Balsara2016.
 */
class WENOAO53 : public Stencil<WENOAO53> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

  // Parameters to control the order reduction (recommended range [0.85,0.95])
  static constexpr double gamma_hi_ = 0.85;
  static constexpr double gamma_lo_ = 0.85;

  // Coefficients for WENOAO53 scheme
  static constexpr double coef_smoothness_1_ = 1.0;
  static constexpr double coef_smoothness_2_ = 13.0 / 3.0;
  static constexpr double coef_smoothness_3_ = 781.0 / 20.0;
  static constexpr double coef_smoothness_4_ = 1421461.0 / 2275.0;

  static constexpr double coef_smoothness_1_1_ = 1.0 / 2.0;
  static constexpr double coef_smoothness_1_2_ = -2.0;
  static constexpr double coef_smoothness_1_3_ = 3.0 / 2.0;
  static constexpr double coef_smoothness_1_4_ = 1.0 / 2.0;
  static constexpr double coef_smoothness_1_5_ = -1.0;
  static constexpr double coef_smoothness_1_6_ = 1.0 / 2.0;

  static constexpr double coef_smoothness_2_1_ = -1.0 / 2.0;
  static constexpr double coef_smoothness_2_2_ = 1.0 / 2.0;
  static constexpr double coef_smoothness_2_3_ = 1.0 / 2.0;
  static constexpr double coef_smoothness_2_4_ = -1.0;
  static constexpr double coef_smoothness_2_5_ = 1.0 / 2.0;

  static constexpr double coef_smoothness_3_1_ = -3.0 / 2.0;
  static constexpr double coef_smoothness_3_2_ = 2.0;
  static constexpr double coef_smoothness_3_3_ = -1.0 / 2.0;
  static constexpr double coef_smoothness_3_4_ = 1.0 / 2.0;
  static constexpr double coef_smoothness_3_5_ = -1.0;
  static constexpr double coef_smoothness_3_6_ = 1.0 / 2.0;

  static constexpr double coef_smoothness_5_01_ = 11.0 / 120.0;
  static constexpr double coef_smoothness_5_02_ = -82.0 / 120.0;
  static constexpr double coef_smoothness_5_03_ = 82.0 / 120.0;
  static constexpr double coef_smoothness_5_04_ = -11.0 / 120.0;
  static constexpr double coef_smoothness_5_05_ = -3.0 / 56.0;
  static constexpr double coef_smoothness_5_06_ = 40.0 / 56.0;
  static constexpr double coef_smoothness_5_07_ = -74.0 / 56.0;
  static constexpr double coef_smoothness_5_08_ = 40.0 / 56.0;
  static constexpr double coef_smoothness_5_09_ = -3.0 / 56.0;
  static constexpr double coef_smoothness_5_10_ = -1.0 / 12.0;
  static constexpr double coef_smoothness_5_11_ = 2.0 / 12.0;
  static constexpr double coef_smoothness_5_12_ = -2.0 / 12.0;
  static constexpr double coef_smoothness_5_13_ = 1.0 / 12.0;
  static constexpr double coef_smoothness_5_14_ = 1.0 / 24.0;
  static constexpr double coef_smoothness_5_15_ = -4.0 / 24.0;
  static constexpr double coef_smoothness_5_16_ = 6.0 / 24.0;
  static constexpr double coef_smoothness_5_17_ = -4.0 / 24.0;
  static constexpr double coef_smoothness_5_18_ = 1.0 / 24.0;

  static constexpr double coef_smoothness_weight_5_1_ = 1.0 / 10.0;
  static constexpr double coef_smoothness_weight_5_2_ = 123.0 / 455.0;

  // Linear weights according to Eq. (3.5)
  static constexpr double linear_weight_r3_1_ =
      (1.0 - gamma_hi_) * (1.0 - gamma_lo_) * 0.5;
  static constexpr double linear_weight_r3_2_ = (1.0 - gamma_hi_) * gamma_lo_;
  static constexpr double linear_weight_r3_3_ = linear_weight_r3_1_;
  static constexpr double linear_weight_r5_3_ = gamma_hi_;
  static constexpr double one_over_linear_weight_r5_3_ =
      1.0 / linear_weight_r5_3_;

  // Legendre polynomials evaluated at cell face x = 1/2
  static constexpr double legendre_1_ = 1.0 / 2.0;
  static constexpr double legendre_2_ = 1.0 / 6.0;
  static constexpr double legendre_3_ = 1.0 / 20.0;
  static constexpr double legendre_4_ = 1.0 / 70.0;

  // Precompiled constant
  static constexpr double one_third_ = 1.0 / 3.0;

  // Number of cells required for upwind and downwind stencils, as well as
  // number of cells downstream of the cell
  static constexpr unsigned int stencil_size_ = 6;
  static constexpr unsigned int downstream_stencil_size_ = 2;

  /**
   * @brief Evaluates the stencil according to WENO-AO(5,3) adaptive-order
   * scheme from Balsara (2016). Also See base class.
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

    // Compute Legendre coefficients u according to Eq. (2.3-2.5), Eq. (2.16)
    // amd smoothness indicators beta according to Eq. (2.6) and Eq. (2.19)
    double const u_r3_11 = coef_smoothness_1_1_ * v1 +
                           coef_smoothness_1_2_ * v2 +
                           coef_smoothness_1_3_ * v3;
    double const u_r3_12 = coef_smoothness_1_4_ * v1 +
                           coef_smoothness_1_5_ * v2 +
                           coef_smoothness_1_6_ * v3;
    double const beta_r3_1 = coef_smoothness_1_ * u_r3_11 * u_r3_11 +
                             coef_smoothness_2_ * u_r3_12 * u_r3_12;

    double const u_r3_21 =
        coef_smoothness_2_1_ * v2 + coef_smoothness_2_2_ * v4;
    double const u_r3_22 = coef_smoothness_2_3_ * v2 +
                           coef_smoothness_2_4_ * v3 +
                           coef_smoothness_2_5_ * v4;
    double const beta_r3_2 = coef_smoothness_1_ * u_r3_21 * u_r3_21 +
                             coef_smoothness_2_ * u_r3_22 * u_r3_22;

    double const u_r3_31 = coef_smoothness_3_1_ * v3 +
                           coef_smoothness_3_2_ * v4 +
                           coef_smoothness_3_3_ * v5;
    double const u_r3_32 = coef_smoothness_3_4_ * v3 +
                           coef_smoothness_3_5_ * v4 +
                           coef_smoothness_3_6_ * v5;
    double const beta_r3_3 = coef_smoothness_1_ * u_r3_31 * u_r3_31 +
                             coef_smoothness_2_ * u_r3_32 * u_r3_32;

    double const u_r5_31 =
        coef_smoothness_5_01_ * v1 + coef_smoothness_5_02_ * v2 +
        coef_smoothness_5_03_ * v4 + coef_smoothness_5_04_ * v5;
    double const u_r5_32 =
        coef_smoothness_5_05_ * v1 + coef_smoothness_5_06_ * v2 +
        coef_smoothness_5_07_ * v3 + coef_smoothness_5_08_ * v4 +
        coef_smoothness_5_09_ * v5;
    double const u_r5_33 =
        coef_smoothness_5_10_ * v1 + coef_smoothness_5_11_ * v2 +
        coef_smoothness_5_12_ * v4 + coef_smoothness_5_13_ * v5;
    double const u_r5_34 =
        coef_smoothness_5_14_ * v1 + coef_smoothness_5_15_ * v2 +
        coef_smoothness_5_16_ * v3 + coef_smoothness_5_17_ * v4 +
        coef_smoothness_5_18_ * v5;

    double const beta_r5_3 =
        coef_smoothness_1_ * (u_r5_31 + coef_smoothness_weight_5_1_ * u_r5_33) *
            (u_r5_31 + coef_smoothness_weight_5_1_ * u_r5_33) +
        coef_smoothness_2_ * (u_r5_32 + coef_smoothness_weight_5_2_ * u_r5_34) *
            (u_r5_32 + coef_smoothness_weight_5_2_ * u_r5_34) +
        coef_smoothness_3_ * u_r5_33 * u_r5_33 +
        coef_smoothness_4_ * u_r5_34 * u_r5_34;

    // Compute normalized weights
    // Eq. (3.6)
    double const tau =
        (Abs(beta_r5_3 - beta_r3_1) + Abs(beta_r5_3 - beta_r3_2) +
         Abs(beta_r5_3 - beta_r3_3)) *
        one_third_;

    // Eq. (3.7a) Note: Balsara et al. suggest an epsilon value of 1e-12 to
    // minimize the influence. We use machine precision instead.
    double const a1 =
        linear_weight_r3_1_ *
        (1.0 + (tau * tau) / ((beta_r3_1 + epsilon_) * (beta_r3_1 + epsilon_)));
    double const a2 =
        linear_weight_r3_2_ *
        (1.0 + (tau * tau) / ((beta_r3_2 + epsilon_) * (beta_r3_2 + epsilon_)));
    double const a3 =
        linear_weight_r3_3_ *
        (1.0 + (tau * tau) / ((beta_r3_3 + epsilon_) * (beta_r3_3 + epsilon_)));
    double const a5 =
        linear_weight_r5_3_ *
        (1.0 + (tau * tau) / ((beta_r5_3 + epsilon_) * (beta_r5_3 + epsilon_)));

    double const one_a_sum = 1.0 / (a1 + a2 + a3 + a5);

    double const w1 = a1 * one_a_sum;
    double const w2 = a2 * one_a_sum;
    double const w3 = a3 * one_a_sum;
    double const w5 = a5 * one_a_sum;

    // Compute coefficients of the Legendre basis polynomial according to Eq.
    // (3.10)
    double const tmp1 = w5 * one_over_linear_weight_r5_3_;

    double const u0 = v3;
    double const u1 =
        tmp1 * (u_r5_31 - linear_weight_r3_1_ * u_r3_11 -
                linear_weight_r3_2_ * u_r3_21 - linear_weight_r3_3_ * u_r3_31) +
        w1 * u_r3_11 + w2 * u_r3_21 + w3 * u_r3_31;
    double const u2 =
        tmp1 * (u_r5_32 - linear_weight_r3_1_ * u_r3_12 -
                linear_weight_r3_2_ * u_r3_22 - linear_weight_r3_3_ * u_r3_32) +
        w1 * u_r3_12 + w2 * u_r3_22 + w3 * u_r3_32;
    double const u3 = tmp1 * u_r5_33;
    double const u4 = tmp1 * u_r5_34;

    // Return value of reconstructed polynomial according to Eq. (3.11) Note:
    // The polynomial constructed in the paper is always evaluated for x = 0.5
    // (cell face).
    return u0 + legendre_1_ * u1 + legendre_2_ * u2 + legendre_3_ * u3 +
           legendre_4_ * u4;
  }

public:
  explicit constexpr WENOAO53() = default;
  ~WENOAO53() = default;
  WENOAO53(WENOAO53 const &) = delete;
  WENOAO53 &operator=(WENOAO53 const &) = delete;
  WENOAO53(WENOAO53 &&) = delete;
  WENOAO53 &operator=(WENOAO53 &&) = delete;
};

#endif // STENCIL_WENOAO53_H
