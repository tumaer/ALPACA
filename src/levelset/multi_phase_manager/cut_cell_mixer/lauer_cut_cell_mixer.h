//===---------------------- lauer_cut_cell_mixer.h ------------------------===//
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
#ifndef LAUER_CUT_CELL_MIXER_H
#define LAUER_CUT_CELL_MIXER_H

#include "two_phase_cut_cell_mixer.h"

/**
 * @brief Performs mixing by a cell-face-aperture based splitting of the mixing
 * fluxes.
 */
class LauerCutCellMixer : public TwoPhaseCutCellMixer<LauerCutCellMixer> {

  friend TwoPhaseCutCellMixer;

  void CalculateMixingContributionsImplementation(
      Node const &node, MaterialName const material,
      std::vector<std::pair<std::vector<std::array<unsigned int, 6>>,
                            std::vector<std::array<double, 2>>>>
          &mixing_contributions) const;

public:
  LauerCutCellMixer() = delete;
  explicit LauerCutCellMixer(HaloManager &halo_manager,
                             MaterialManager const &material_manager);
  ~LauerCutCellMixer() = default;
  LauerCutCellMixer(LauerCutCellMixer const &) = delete;
  LauerCutCellMixer &operator=(LauerCutCellMixer const &) = delete;
  LauerCutCellMixer(LauerCutCellMixer &&) = delete;
  LauerCutCellMixer &operator=(LauerCutCellMixer &&) = delete;
};

#endif // LAUER_CUT_CELL_MIXER_H
