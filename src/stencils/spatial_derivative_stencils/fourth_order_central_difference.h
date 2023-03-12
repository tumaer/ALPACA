//===---------------- fourth_order_central_difference.h -------------------===//
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
#ifndef FOURTH_ORDER_CENTRAL_DIFFERENCE_H
#define FOURTH_ORDER_CENTRAL_DIFFERENCE_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialDerivativeStencil class to evaluate the
 * stencil with a 4th order central differencing scheme. See also base class.
 */
class FourthOrderCentralDifference
    : public Stencil<FourthOrderCentralDifference> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Derivative;

  static constexpr unsigned int stencil_size_ = 5;
  static constexpr unsigned int downstream_stencil_size_ = 2;

  /**
   * @brief Evaluates the stencil according to a fourth order central scheme.
   * Also See base class.
   * @note Hotpath function.
   */
  constexpr double
  ApplyImplementation(std::array<double, stencil_size_> const &array,
                      std::array<int const, 2> const,
                      double const cell_size) const {
    double const denominator = cell_size * 12.0;
    double const result = -array[downstream_stencil_size_ + 2] +
                          8.0 * array[downstream_stencil_size_ + 1] -
                          8.0 * array[downstream_stencil_size_ - 1] +
                          array[downstream_stencil_size_ - 2];
    return result / denominator;
  }

public:
  explicit constexpr FourthOrderCentralDifference() = default;
  ~FourthOrderCentralDifference() = default;
};

#endif // FOURTH_ORDER_CENTRAL_DIFFERENCE_H
