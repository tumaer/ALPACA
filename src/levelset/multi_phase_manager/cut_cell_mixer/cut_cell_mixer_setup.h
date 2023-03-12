//===---------------------- cut_cell_mixer_setup.h ------------------------===//
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
#ifndef CUT_CELL_MIXER_SETUP_H
#define CUT_CELL_MIXER_SETUP_H

#include "aperture_cut_cell_mixer.h"
#include "lauer_cut_cell_mixer.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a CutCellMixer type based on a specified constexpr.
 */
namespace CutCellMixerSetup {

/**
 * @brief Function returning the typedef of a CutCellMixer based on a constexpr
 * template.
 *
 * @tparam CutCellMixers The constexpr template parameter to specify the exact
 * CutCellMixer type.
 */
template <CutCellMixers> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<CutCellMixers::ApertureBased> {
  typedef ApertureCutCellMixer type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<CutCellMixers::Lauer> {
  typedef LauerCutCellMixer type;
};

} // namespace CutCellMixerSetup

#endif // CUT_CELL_MIXER_SETUP_H
