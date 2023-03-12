//===------------------------ cut_cell_mixer.h ----------------------------===//
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
#ifndef CUT_CELL_MIXER_H
#define CUT_CELL_MIXER_H

#include "halo_manager.h"
#include "levelset/geometry/geometry_calculator_setup.h"
#include "materials/material_manager.h"
#include "user_specifications/numerical_setup.h"

using GeometryCalculatorConcretization =
    GeometryCalculatorSetup::Concretize<geometry_calculator>::type;

/**
 * @brief The CutCellMixer class mixes small cut-cells with its neighbors.
 * @tparam DerivedCutCellMixer Typename as template parameter due to CRTP.
 */
template <typename DerivedCutCellMixer> class CutCellMixer {

protected:
  const GeometryCalculatorConcretization geometry_calculator_;
  HaloManager &halo_manager_;
  MaterialManager const &material_manager_;

  /**
   * @brief Default constructor of the CutCellMixer class.
   * @param halo_manager Instance to a HaloManager which provides MPI-related
   * methods.
   */
  explicit CutCellMixer(HaloManager &halo_manager,
                        MaterialManager const &material_manager)
      : geometry_calculator_(), halo_manager_(halo_manager),
        material_manager_(material_manager) {
    // Empty Constructor, besides initializer list.
  }

public:
  explicit CutCellMixer() = default;
  ~CutCellMixer() = default;
  CutCellMixer(CutCellMixer const &) = delete;
  CutCellMixer &operator=(CutCellMixer const &) = delete;
  CutCellMixer(CutCellMixer &&) = delete;
  CutCellMixer &operator=(CutCellMixer &&) = delete;
  /**
   * @brief Provides functionality for a cut-cell mixing procedure.
   * @param node The node for which mixing has to be performed.
   */
  void Mix(Node &node) const {
    static_cast<DerivedCutCellMixer const &>(*this).MixImplementation(node);
  }
};

#endif // CUT_CELL_MIXER_H
