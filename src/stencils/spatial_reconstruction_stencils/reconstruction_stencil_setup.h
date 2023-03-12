//===----------------- reconstruction_stencil_setup.h ---------------------===//
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
#ifndef RECONSTRUCTION_STENCIL_SETUP_H
#define RECONSTRUCTION_STENCIL_SETUP_H

#include "first_order.h"
#include "fourth_order_central.h"
#include "teno5.h"
#include "user_specifications/stencil_setup.h"
#include "weno-ao53.h"
#include "weno-cu6.h"
#include "weno3.h"
#include "weno5.h"
#include "weno5_hm.h"
#include "weno5_is.h"
#include "weno5_nu6p.h"
#include "weno5_z.h"
#include "weno7.h"
#include "weno9.h"
#include "weno_f3p.h"

/**
 * @brief A namespace to get a ReconstructionStencil type based on a specified
 * constexpr.
 */
namespace ReconstructionStencilSetup {

/**
 * @brief Function returning the typedef of a ReconstructionStencil based on a
 * constexpr template.
 *
 * @tparam ReconstructionStencils The constexpr template parameter to specify
 * the exact ReconstructionStencil type.
 */
template <ReconstructionStencils> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::FirstOrder> {
  typedef FirstOrder type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENO3> {
  typedef WENO3 type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENOF3P> {
  typedef WENOF3P type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::FourthOrderCentral> {
  typedef FourthOrderCentral type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENO5> {
  typedef WENO5 type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENO5IS> {
  typedef WENO5IS type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENO5Z> {
  typedef WENO5Z type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENOAO53> {
  typedef WENOAO53 type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENO5NU6P> {
  typedef WENO5NU6P type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::TENO5> {
  typedef TENO5 type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENOCU6> {
  typedef WENOCU6 type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENO7> {
  typedef WENO7 type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENO9> {
  typedef WENO9 type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ReconstructionStencils::WENO5HM> {
  typedef WENO5HM type;
};

} // namespace ReconstructionStencilSetup

#endif // RECONSTRUCTION_STENCIL_SETUP_H
