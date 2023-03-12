//===----------------------- parameter_manager.cpp ------------------------===//
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
#include "parameter/parameter_manager.h"

#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief Standard constructor to create the ParameterManager object.
 * @param material_manager Instance for handling models of material proeprties
 * (such as viscosity).
 * @param halo_manager Instance to provide halo updates.
 */
ParameterManager::ParameterManager(MaterialManager const &material_manager,
                                   HaloManager &halo_manager)
    : material_manager_(material_manager),
      ghost_fluid_extender_(material_manager, halo_manager) {
  /** Empty besides initializer list */
}

/**
 * @brief Updates all parameters of a single node.
 * @param node The node under consideration (indirect return).
 */
void ParameterManager::UpdateParameters(Node &node) const {

  // Obtain the cell size for the given node
  double const cell_size = node.GetCellSize();

  // Update all parameters acting on a single material
  for (auto &phase : node.GetPhases()) {

    // Obtain the sign of the material and the material
    std::int8_t const material_sign =
        MaterialSignCapsule::SignOfMaterial(phase.first);
    Material const &material = material_manager_.GetMaterial(phase.first);

    // Start block to update the material property models
    // Shear viscosity
    if constexpr (CC::ViscosityIsActive() && CC::ShearViscosityModelActive()) {
      // Call appropriate function depending on presence of an interface of this
      // node
      if (node.HasLevelset()) {
        material.GetShearViscosityModel().UpdateParameter(
            phase.second, cell_size,
            node.GetInterfaceTags<
                InterfaceDescriptionBufferType::Reinitialized>(),
            material_sign);
      } else {
        material.GetShearViscosityModel().UpdateParameter(phase.second,
                                                          cell_size);
      }
    }

    // Thermal conductivity
    if constexpr (CC::HeatConductionActive() &&
                  CC::ThermalConductivityModelActive()) {
      // Call appropriate function depending on presence of an interface of this
      // node
      if (node.HasLevelset()) {
        material.GetThermalConductivityModel().UpdateParameter(
            phase.second, cell_size,
            node.GetInterfaceTags<
                InterfaceDescriptionBufferType::Reinitialized>(),
            material_sign);
      } else {
        material.GetThermalConductivityModel().UpdateParameter(phase.second,
                                                               cell_size);
      }
    }
  }
}

/**
 * @brief Extends for all nodes the parameters into the extension band.
 * @param nodes Nodes on which the extension should be performed.
 */
void ParameterManager::ExtendParameters(
    std::vector<std::reference_wrapper<Node>> const &nodes) const {
  ghost_fluid_extender_.Extend(nodes);
}
