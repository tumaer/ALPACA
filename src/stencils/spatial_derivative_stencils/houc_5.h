//===---------------------------- houc_5.h --------------------------------===//
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
#ifndef HOUC_5_H
#define HOUC_5_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialDerivativeStencil class to evaluate the
 * stencil with a 5th order HOUC stencil of \cite Nourgaliev2007.
 */
class HOUC5 : public Stencil<HOUC5> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Derivative;

  static constexpr unsigned int stencil_size_ = 7;
  static constexpr unsigned int downstream_stencil_size_ = 3;

  static constexpr double one_sixtieth = 1.0 / 60.0;

  static constexpr double coefficient_0_ = -2.0;
  static constexpr double coefficient_1_ = +15.0;
  static constexpr double coefficient_2_ = -60.0;
  static constexpr double coefficient_3_ = +20.0;
  static constexpr double coefficient_4_ = +30.0;
  static constexpr double coefficient_5_ = -3.0;

  /**
   * @brief Implements the 5th order HOUC stencil. Also See base class.
   * @note Hotpath function.
   */
  constexpr double
  ApplyImplementation(std::array<double, stencil_size_> const &array,
                      std::array<int const, 2> const evaluation_properties,
                      double const) const {
    // Assign values to v_i to make it easier to read
    double const v0 =
        array[downstream_stencil_size_ - 3 * evaluation_properties[1]];
    double const v1 =
        array[downstream_stencil_size_ - 2 * evaluation_properties[1]];
    double const v2 =
        array[downstream_stencil_size_ - 1 * evaluation_properties[1]];
    double const v3 = array[downstream_stencil_size_];
    double const v4 =
        array[downstream_stencil_size_ + 1 * evaluation_properties[1]];
    double const v5 =
        array[downstream_stencil_size_ + 2 * evaluation_properties[1]];

    double const result = coefficient_0_ * v0 + coefficient_1_ * v1 +
                          coefficient_2_ * v2 + coefficient_3_ * v3 +
                          coefficient_4_ * v4 + coefficient_5_ * v5;

    return result * one_sixtieth * evaluation_properties[1];
  }

public:
  explicit constexpr HOUC5() = default;
  ~HOUC5() = default;
  HOUC5(HOUC5 const &) = delete;
  HOUC5 &operator=(HOUC5 const &) = delete;
  HOUC5(HOUC5 &&) = delete;
  HOUC5 &operator=(HOUC5 &&) = delete;
};

#endif // HOUC_5_H
