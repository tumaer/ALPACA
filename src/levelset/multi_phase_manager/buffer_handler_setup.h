//===---------------------- buffer_handler_setup.h ------------------------===//
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
#ifndef BUFFER_HANDLER_SETUP_H
#define BUFFER_HANDLER_SETUP_H

#include "two_phase_buffer_handler.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a BufferHandler type based on a specified
 * constexpr.
 */
namespace BufferHandlerSetup {

/**
 * @brief Function returning the typedef of a BufferHandler based on a constexpr
 * template.
 *
 * @tparam BufferHandlers The constexpr template parameter to specify the exact
 * BufferHandler type.
 */
template <BufferHandlers> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<BufferHandlers::TwoPhase> {
  typedef TwoPhaseBufferHandler type;
};

} // namespace BufferHandlerSetup

#endif // BUFFER_HANDLER_SETUP_H
