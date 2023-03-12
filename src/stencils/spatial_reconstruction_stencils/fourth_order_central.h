//===---------------------- fourth_order_central.h ------------------------===//
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
#ifndef FOURTH_ORDER_CENTRAL_H
#define FOURTH_ORDER_CENTRAL_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to evaluate
 * the stencil with a fourth-order central scheme.
 */
class FourthOrderCentral : public Stencil<FourthOrderCentral> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

  static constexpr double one_sixteenth_ = 1.0 / 16.0;

  // Number of cells required for upwind and downwind stencils, as well as
  // number of cells downstream of the cell
  static constexpr unsigned int stencil_size_ = 4;
  static constexpr unsigned int downstream_stencil_size_ = 1;

  /**
   * @brief Evaluates the stencil according to a fourth order central scheme.
   * Also See base class.
   * @note Hotpath function.
   */
  constexpr double
  ApplyImplementation(std::array<double, stencil_size_> const &array,
                      std::array<int const, 2> const, double const) const {
    double const result = 9.0 * (array[downstream_stencil_size_ - 0] +
                                 array[downstream_stencil_size_ + 1]) -
                          1.0 * (array[downstream_stencil_size_ - 1] +
                                 array[downstream_stencil_size_ + 2]);
    return result * one_sixteenth_;
  }

public:
  explicit constexpr FourthOrderCentral() = default;
  ~FourthOrderCentral() = default;
};

#endif // STENCIL_FOURTH_ORDER_CENTRAL_H
