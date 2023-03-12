//===----------------------- central_difference.h -------------------------===//
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
#ifndef CENTRAL_DIFFERENCE_H
#define CENTRAL_DIFFERENCE_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialDerivativeStencil class to evaluate the
 * stencil with a 2nd order central differencing scheme. See also base class.
 */
class CentralDifference : public Stencil<CentralDifference> {

  friend Stencil;

  static constexpr StencilType stencil_type_ = StencilType::Derivative;

  static constexpr unsigned int stencil_size_ = 3;
  static constexpr unsigned int downstream_stencil_size_ = 1;

  /**
   * @brief Evaluates the stencil according to a fourth order central scheme.
   * Also See base class.
   * @note Hotpath function.
   */
  constexpr double
  ApplyImplementation(std::array<double, stencil_size_> const &array,
                      std::array<int const, 2> const,
                      double const cell_size) const {
    return 0.5 *
           (array[downstream_stencil_size_ + 1] -
            array[downstream_stencil_size_ - 1]) /
           cell_size;
  }

public:
  explicit constexpr CentralDifference() = default;
  ~CentralDifference() = default;
};

#endif // CENTRAL_DIFFERENCE_H
