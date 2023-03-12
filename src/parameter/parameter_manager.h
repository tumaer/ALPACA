//===----------------------- parameter_manager.h --------------------------===//
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
#ifndef PARAMETER_MANAGER_H
#define PARAMETER_MANAGER_H

#include "levelset/multi_phase_manager/ghost_fluid_extender/ghost_fluid_extender_setup.h"
#include "materials/material_manager.h"
#include "parameter/material_parameter_model.h"
#include "topology/node.h"

using GhostFluidExtenderConcretizationParameter =
    GhostFluidExtenderSetup::Concretize<extender>::type_parameters;

/**
 * @brief The ParameterManager handles all updates regarding the parameter
 * buffers that lie on a single material block. In general, it is not restricted
 * to models, such as viscosity or thermal conductivity models, that act on
 * material properties. It can handle all models that act on a specific
 * parameter buffer.
 */
class ParameterManager {
  // Instance for handling parameter calculation of material data (e.g.,
  // viscosity)
  MaterialManager const &material_manager_;
  // Instance to extend parameters
  GhostFluidExtenderConcretizationParameter const ghost_fluid_extender_;

public:
  ParameterManager() = delete;
  explicit ParameterManager(MaterialManager const &material_manager,
                            HaloManager &halo_manager);
  ~ParameterManager() = default;
  ParameterManager(ParameterManager const &) = delete;
  ParameterManager &operator=(ParameterManager const &) = delete;
  ParameterManager(ParameterManager &&) = delete;
  ParameterManager &operator=(ParameterManager &&) = delete;

  // Update functions for all parameters
  void UpdateParameters(Node &node) const;
  // Extension function for all parameters
  void ExtendParameters(
      std::vector<std::reference_wrapper<Node>> const &nodes) const;
};

#endif // PARAMETER_MANAGER_H
