//===-------------------------- first_order.h -----------------------------===//
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
#ifndef FIRST_ORDER_H
#define FIRST_ORDER_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to evaluate
 * the stencil with a first-order scheme.
 */
class FirstOrder : public Stencil<FirstOrder> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

  // Number of cells required for upwind and downwind stencils, as well as
  // number of cells downstream of the cell
  static constexpr unsigned int stencil_size_ = 2;
  static constexpr unsigned int downstream_stencil_size_ = 0;

  /**
   * @brief Evaluates the stencil according to a first order scheme. Also See
   * base class.
   * @note Hotpath function.
   */
  constexpr double
  ApplyImplementation(std::array<double, stencil_size_> const &array,
                      std::array<int const, 2> const evaluation_properties,
                      double const) const {
    // Return left/right value
    return array[evaluation_properties[0]];
  }

public:
  explicit constexpr FirstOrder() = default;
  ~FirstOrder() = default;
  FirstOrder(FirstOrder const &) = delete;
  FirstOrder &operator=(FirstOrder const &) = delete;
  FirstOrder(FirstOrder &&) = delete;
  FirstOrder &operator=(FirstOrder &&) = delete;
};

#endif // STENCIL_FIRST_ORDER_H
