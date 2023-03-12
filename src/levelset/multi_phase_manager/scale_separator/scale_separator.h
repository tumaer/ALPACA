//===----------------------- scale_separator.h ----------------------------===//
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
#ifndef SCALE_SEPARATOR_H
#define SCALE_SEPARATOR_H

#include "halo_manager.h"
#include "levelset/geometry/geometry_calculator_marching_cubes.h"
#include "materials/material_manager.h"
#include "user_specifications/numerical_setup.h"

template <typename DerivedScaleSeparator> class ScaleSeparator {

  friend DerivedScaleSeparator;

protected:
  MaterialManager const &material_manager_;
  HaloManager
      &halo_manager_; // TODO-19 NH Think about making it const (rats tail)

  ScaleSeparator() = delete;

  /**
   * @brief Default constructor of the ScaleSeparator class.
   * @param material_manager Instance of a MaterialManager, which already has
   * been initialized according to the user input.
   * @param halo_manager Instance to a HaloManager which provides MPI-related
   * methods.
   */
  explicit ScaleSeparator(MaterialManager const &material_manager,
                          HaloManager &halo_manager)
      : material_manager_(material_manager), halo_manager_(halo_manager) {
    // Empty Constructor, besides initializer list.
  }

public:
  virtual ~ScaleSeparator() = default;

  /**
   * @brief Performs a scale separation procedure.
   * @param nodes The nodes for which scale separation should be done.
   * @param buffer_type The level-set buffer type for which scale separation is
   * done.
   */
  void SeparateScales(std::vector<std::reference_wrapper<Node>> const &nodes,
                      InterfaceBlockBufferType const buffer_type) const {
    static_cast<DerivedScaleSeparator const &>(*this)
        .SeparateScalesImplementation(nodes, buffer_type);
  }
};

#endif // SCALE_SEPARATOR_H
