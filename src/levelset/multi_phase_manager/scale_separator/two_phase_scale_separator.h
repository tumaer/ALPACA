//===-------------------- two_phase_scale_separator.h ---------------------===//
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
#ifndef TWO_PHASE_SCALE_SEPARATOR_H
#define TWO_PHASE_SCALE_SEPARATOR_H

#include "scale_separator.h"

/**
 * @brief Class providing the functionality to identify and handle scale
 * separation for levelset fields for two different phases.
 */
class TwoPhaseScaleSeparator : public ScaleSeparator<TwoPhaseScaleSeparator> {

  friend ScaleSeparator;

private:
  void SeparateScalesImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes,
      InterfaceBlockBufferType const buffer_type) const;

public:
  TwoPhaseScaleSeparator() = delete;
  explicit TwoPhaseScaleSeparator(MaterialManager const &material_manager,
                                  HaloManager &halo_manager);
  ~TwoPhaseScaleSeparator() = default;
};

#endif // TWO_PHASE_SCALE_SEPARATOR_H
