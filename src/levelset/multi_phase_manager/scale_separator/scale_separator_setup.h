//===-------------------- scale_separator_setup.h -------------------------===//
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
#ifndef SCALE_SEPARATOR_SETUP_H
#define SCALE_SEPARATOR_SETUP_H

#include "two_phase_scale_separator.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a ScaleSeparator type based on a specified
 * constexpr.
 */
namespace ScaleSeparatorSetup {

/**
 * @brief Function returning the typedef of a ScaleSeparator based on a
 * constexpr template.
 *
 * @tparam ScaleSeparators The constexpr template parameter to specify the exact
 * ScaleSeparator type.
 */
template <ScaleSeparators> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<ScaleSeparators::TwoPhase> {
  typedef TwoPhaseScaleSeparator type;
};

} // namespace ScaleSeparatorSetup

#endif // SCALE_SEPARATOR_SETUP_H
