//===-------------------------- weno5_nu6p.h ------------------------------===//
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
#ifndef WENO5NU6P_H
#define WENO5NU6P_H

#include "stencils/stencil.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute
 * fluxes according to WENO5-nu6+ \cite Gande2020
 */
class WENO5NU6P : public Stencil<WENO5NU6P> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

  // Coefficients for WENO5Z scheme
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

  static constexpr double coef_nu6_ = 1.0 / 144.0;

  // Number of cells required for upwind and downwind stencils, as well as
  // number of cells downstream of the cell
  static constexpr unsigned int stencil_size_ = 6;
  static constexpr unsigned int downstream_stencil_size_ = 2;

  /**
   * @brief Evaluates the stencil according to a WENO5-NU6P scheme. Also See
   * base class.
   * @note Hotpath function.
   */
  constexpr double
  ApplyImplementation(std::array<double, stencil_size_> const &array,
                      std::array<int const, 2> const evaluation_properties,
                      double const cell_size) const {
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
    double const s11 = (v1 + 3.0 * v3) - 4.0 * v2;
    double const s12 = (v1 + v3) - 2.0 * v2;

    double nu0 = 0.25 * (s11 * s11) + (s12 * s12);

    double const s21 = v2 - v4;
    double const s22 = (v2 + v4) - 2.0 * v3;
    double const nu1 = 0.25 * (s21 * s21) + (s22 * s22);

    double const s31 = (v5 + 3.0 * v3) - 4.0 * v4;
    double const s32 = (v5 + v3) - 2.0 * v4;
    double const nu2 = 0.25 * (s31 * s31) + (s32 * s32);

    double const s41 = v1 - 8.0 * v2 + 8.0 * v3 - v5;
    double const s42 = v1 - 16.0 * v2 + 30.0 * v3 - 16.0 * v4 + v5;
    const double nu5 = coef_nu6_ * (s41 * s41) + (s42 * s42);

    // Compute weights Note: Borges et al. suggest an epsilon value of 1e-40 to
    // minimize the influence. We use machine precision instead.
    double const tau = std::abs(nu5 - (4.0 * nu1 + nu0 + nu2) / 6.0);
    // Note: Gande et al. suggest lambda = cell_size ^ 101/102. We use a cheap
    // unity exponent instead.
    double const lambda = cell_size;

    double const tmp = lambda / (tau + epsilon_);
    double const a0 =
        coef_weights_1_ *
        (1.0 + (tau + epsilon_) / (nu0 + epsilon_) + tmp * (nu0 + epsilon_));
    double const a1 =
        coef_weights_2_ *
        (1.0 + (tau + epsilon_) / (nu1 + epsilon_) + tmp * (nu1 + epsilon_));
    double const a2 =
        coef_weights_3_ *
        (1.0 + (tau + epsilon_) / (nu2 + epsilon_) + tmp * (nu2 + epsilon_));

    double const one_a_sum = 1.0 / (a0 + a1 + a2);

    double const w1 = a0 * one_a_sum;
    double const w2 = a1 * one_a_sum;
    double const w3 = a2 * one_a_sum;

    // Return weighted average
    return w1 * (coef_stencils_1_ * v1 + coef_stencils_2_ * v2 +
                 coef_stencils_3_ * v3) +
           w2 * (coef_stencils_4_ * v2 + coef_stencils_5_ * v3 +
                 coef_stencils_6_ * v4) +
           w3 * (coef_stencils_7_ * v3 + coef_stencils_8_ * v4 +
                 coef_stencils_9_ * v5);
  }

public:
  explicit constexpr WENO5NU6P() = default;
  ~WENO5NU6P() = default;
  WENO5NU6P(WENO5NU6P const &) = delete;
  WENO5NU6P &operator=(WENO5NU6P const &) = delete;
  WENO5NU6P(WENO5NU6P &&) = delete;
  WENO5NU6P &operator=(WENO5NU6P &&) = delete;
};

#endif // STENCIL_WENO5NU6P_H
