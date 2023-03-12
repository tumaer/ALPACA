//===--------------------------- weno5_is.h -------------------------------===//
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
#ifndef WENO5IS_H
#define WENO5IS_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute
 * fluxes according to the incremental-stencil WENO \cite Wang2018
 */
class WENO5IS : public Stencil<WENO5IS> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

  // Coefficients for WENO5IS scheme
  static constexpr double d0_ = 15.0 / 32.0;
  static constexpr double d1_ = 5.0 / 32.0;
  static constexpr double d2_ = 5.0 / 16.0;
  static constexpr double d3_ = 1.0 / 16.0;

  static constexpr double c1_ = 13.0 / 12.0;

  // Small values to avoid division by 0, but also to adjust dissipation.
  static constexpr double epsilon_weno5_ = 1.0e-6;

  // Number of cells required for upwind and downwind stencils, as well as
  // number of cells downstream of the cell
  static constexpr unsigned int stencil_size_ = 6;
  static constexpr unsigned int downstream_stencil_size_ = 2;

  /**
   * @brief Evaluates the stencil according to a WENO5IS scheme. Also See base
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

    // Compute smoothness indicators beta^_i
    double const tmp1 = v3 - v4;
    double const beta0 = tmp1 * tmp1;

    double const tmp2 = v2 - v3;
    double const beta1 = tmp2 * tmp2;

    double const tmp3 = v2 - 2.0 * v3 + v4;
    double const tmp4 = v2 - v4;
    double const beta01 = c1_ * tmp3 * tmp3 + 0.25 * tmp4 * tmp4;

    double const tmp5 = v3 - 2.0 * v4 + v5;
    double const tmp6 = 3.0 * v3 - 4.0 * v4 + v5;
    double const beta2 = c1_ * tmp5 * tmp5 + 0.25 * tmp6 * tmp6;

    double const tmp7 = v1 - 2.0 * v2 + v3;
    double const tmp8 = v1 - 4.0 * v2 + 3.0 * v3;
    double const beta3 = c1_ * tmp7 * tmp7 + 0.25 * tmp8 * tmp8;

    double const tmp9 = v5 - 4.0 * v4 + 6.0 * v3 - 4.0 * v2 + v1;
    double const tmp10 = v5 - 2.0 * v4 + 2.0 * v2 - v1;
    double const tau5 = c1_ * tmp9 * tmp9 + 0.25 * tmp10 * tmp10;

    // Compute weights
    double const a0 = d0_ * (1.0 + (tau5 / (beta0 + epsilon_)) *
                                       (tau5 / (beta01 + epsilon_)));
    double const a1 = d1_ * (1.0 + (tau5 / (beta1 + epsilon_)) *
                                       (tau5 / (beta01 + epsilon_)));
    double const a2 = d2_ * (1.0 + tau5 / (beta2 + epsilon_));
    double const a3 = d3_ * (1.0 + tau5 / (beta3 + epsilon_));
    double const one_a_sum = 1.0 / (a0 + a1 + a2 + a3);

    double const w0 = a0 * one_a_sum;
    double const w1 = a1 * one_a_sum;
    double const w2 = a2 * one_a_sum;
    double const w3 = a3 * one_a_sum;

    double const u0 = 0.5 * (v3 + v4);
    double const u1 = 0.5 * (3.0 * v3 - v2);
    double const u2 = 0.125 * (3.0 * v3 + 6.0 * v4 - v5);
    double const u3 = 0.125 * (3.0 * v1 - 10.0 * v2 + 15.0 * v3);

    // Return weighted average
    return w0 * u0 + w1 * u1 + w2 * u2 + w3 * u3;
  }

public:
  explicit constexpr WENO5IS() = default;
  ~WENO5IS() = default;
  WENO5IS(WENO5IS const &) = delete;
  WENO5IS &operator=(WENO5IS const &) = delete;
  WENO5IS(WENO5IS &&) = delete;
  WENO5IS &operator=(WENO5IS &&) = delete;
};

#endif // STENCIL_WENO5_IS_H
