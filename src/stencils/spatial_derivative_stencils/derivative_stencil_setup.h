//===------------------ derivative_stencil_setup.h ------------------------===//
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
#ifndef DERIVATIVE_STENCIL_SETUP_H
#define DERIVATIVE_STENCIL_SETUP_H

#include "central_difference.h"
#include "fourth_order_cell_face.h"
#include "fourth_order_central_difference.h"
#include "houc_5.h"
#include "user_specifications/stencil_setup.h"

/**
 * @brief A namespace to get a DerivativeStencil type based on a specified
 * constexpr.
 */
namespace DerivativeStencilSetup {

/**
 * @brief Function returning the typedef of a DerivativeStencil based on a
 * constexpr template.
 *
 * @tparam DerivativeStencils The constexpr template parameter to specify the
 * exact DerivativeStencil type.
 */
template <DerivativeStencils> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<DerivativeStencils::CentralDifference> {
  typedef CentralDifference type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<DerivativeStencils::FourthOrderCellFace> {
  typedef FourthOrderCellFace type;
};
/**
 * @brief See generic implementation.
 */
template <>
struct Concretize<DerivativeStencils::FourthOrderCentralDifference> {
  typedef FourthOrderCentralDifference type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<DerivativeStencils::HOUC5> {
  typedef HOUC5 type;
};

} // namespace DerivativeStencilSetup

#endif // DERIVATIVE_STENCIL_SETUP_H
