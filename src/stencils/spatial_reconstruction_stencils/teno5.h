//===---------------------------- teno5.h ---------------------------------===//
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
#ifndef TENO5_H
#define TENO5_H

#include "stencils/stencil.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute
 * fluxes according to \cite Fu2016a.
 */
class TENO5 : public Stencil<TENO5> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

  // Ideal weights and threshold for stencil smoothness of TENO5
  static constexpr double d1_ = 0.6;
  static constexpr double d2_ = 0.3;
  static constexpr double d3_ = 0.1;
  static constexpr double CT_ = 1.0e-5;

  // Coefficients for TENO5 scheme
  static constexpr double coef_smoothness_1_ = 13.0 / 12.0;
  static constexpr double coef_smoothness_2_ = 0.25;

  static constexpr double coef_smoothness_11_ = 1.0;
  static constexpr double coef_smoothness_12_ = -2.0;
  static constexpr double coef_smoothness_13_ = 1.0;
  static constexpr double coef_smoothness_14_ = 1.0;
  static constexpr double coef_smoothness_15_ = -1.0;

  static constexpr double coef_smoothness_21_ = 1.0;
  static constexpr double coef_smoothness_22_ = -2.0;
  static constexpr double coef_smoothness_23_ = 1.0;
  static constexpr double coef_smoothness_24_ = 3.0;
  static constexpr double coef_smoothness_25_ = -4.0;
  static constexpr double coef_smoothness_26_ = 1.0;

  static constexpr double coef_smoothness_31_ = 1.0;
  static constexpr double coef_smoothness_32_ = -2.0;
  static constexpr double coef_smoothness_33_ = 1.0;
  static constexpr double coef_smoothness_34_ = 1.0;
  static constexpr double coef_smoothness_35_ = -4.0;
  static constexpr double coef_smoothness_36_ = 3.0;

  static constexpr double coef_stencils_1_ = -1.0;
  static constexpr double coef_stencils_2_ = 5.0;
  static constexpr double coef_stencils_3_ = 2.0;
  static constexpr double coef_stencils_4_ = 2.0;
  static constexpr double coef_stencils_5_ = 5.0;
  static constexpr double coef_stencils_6_ = -1.0;
  static constexpr double coef_stencils_7_ = 2.0;
  static constexpr double coef_stencils_8_ = -7.0;
  static constexpr double coef_stencils_9_ = 11.0;

  static constexpr double multiplyer_stencils_ = 1.0 / 6.0;

  // Number of cells required for upwind and downwind stencils, as well as
  // number of cells downstream of the cell
  static constexpr unsigned int stencil_size_ = 6;
  static constexpr unsigned int downstream_stencil_size_ = 2;

  /**
   * @brief Evaluates the stencil according to a TENO scheme of fifth order.
   * Also See base class.
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

    // Compute smoothness indicators si
    double const s11 = coef_smoothness_11_ * v2 + coef_smoothness_12_ * v3 +
                       coef_smoothness_13_ * v4;
    double const s12 = coef_smoothness_14_ * v2 + coef_smoothness_15_ * v4;

    double const s1 =
        coef_smoothness_1_ * s11 * s11 + coef_smoothness_2_ * s12 * s12;

    double const s21 = coef_smoothness_21_ * v3 + coef_smoothness_22_ * v4 +
                       coef_smoothness_23_ * v5;
    double const s22 = coef_smoothness_24_ * v3 + coef_smoothness_25_ * v4 +
                       coef_smoothness_26_ * v5;

    double const s2 =
        coef_smoothness_1_ * s21 * s21 + coef_smoothness_2_ * s22 * s22;

    double const s31 = coef_smoothness_31_ * v1 + coef_smoothness_32_ * v2 +
                       coef_smoothness_33_ * v3;
    double const s32 = coef_smoothness_34_ * v1 + coef_smoothness_35_ * v2 +
                       coef_smoothness_36_ * v3;

    double const s3 =
        coef_smoothness_1_ * s31 * s31 + coef_smoothness_2_ * s32 * s32;

    double const tau5 = Abs(s3 - s2);

    double a1 = 1.0 + tau5 / (s1 + epsilon_);
    double a2 = 1.0 + tau5 / (s2 + epsilon_);
    double a3 = 1.0 + tau5 / (s3 + epsilon_);

    // NF Calculate a^6 without using std::pow to improve performance
    // drastically
    a1 *= (a1 * a1);
    a2 *= (a2 * a2);
    a3 *= (a3 * a3);
    a1 *= a1;
    a2 *= a2;
    a3 *= a3;

    double const one_a_sum = 1.0 / (a1 + a2 + a3);

    double const b1 = a1 * one_a_sum < CT_ ? 0.0 : 1.0;
    double const b2 = a2 * one_a_sum < CT_ ? 0.0 : 1.0;
    double const b3 = a3 * one_a_sum < CT_ ? 0.0 : 1.0;

    double const Variation1 =
        coef_stencils_1_ * v2 + coef_stencils_2_ * v3 + coef_stencils_3_ * v4;
    double const Variation2 =
        coef_stencils_4_ * v3 + coef_stencils_5_ * v4 + coef_stencils_6_ * v5;
    double const Variation3 =
        coef_stencils_7_ * v1 + coef_stencils_8_ * v2 + coef_stencils_9_ * v3;

    double const w1 = d1_ * b1;
    double const w2 = d2_ * b2;
    double const w3 = d3_ * b3;

    double const one_w_sum = 1.0 / (w1 + w2 + w3);

    double const w1_normalized = w1 * one_w_sum;
    double const w2_normalized = w2 * one_w_sum;
    double const w3_normalized = w3 * one_w_sum;

    return (w1_normalized * Variation1 + w2_normalized * Variation2 +
            w3_normalized * Variation3) *
           multiplyer_stencils_;
  }

public:
  explicit constexpr TENO5() = default;
  ~TENO5() = default;
  TENO5(TENO5 const &) = delete;
  TENO5 &operator=(TENO5 const &) = delete;
  TENO5(TENO5 &&) = delete;
  TENO5 &operator=(TENO5 &&) = delete;
};

#endif // STENCIL_TENO5_H
