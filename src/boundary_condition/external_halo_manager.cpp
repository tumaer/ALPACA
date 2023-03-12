//===------------------- external_halo_manager.cpp ------------------------===//
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
#include "boundary_condition/external_halo_manager.h"

/**
 * @brief Constructs the external halo manager with the given input data.
 * @param material_boundary_conditions Array with already initialized material
 * boundary conditions (in total six, one for each side in 3D) (ownership
 * transfer takes place).
 * @param levelset_boundary_conditions Array with already initialized levelset
 * boundary conditions (in total six, one for each side in 3D) (ownership
 * transfer takes place).
 * @note In case of 1D and 2D simulations the array still has its maximum size,
 * but the not used conditions are nullptr.
 */
ExternalHaloManager::ExternalHaloManager(
    std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6>
        material_boundary_conditions,
    std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6>
        levelset_boundary_conditions)
    : material_boundary_conditions_(std::move(material_boundary_conditions)),
      levelset_boundary_conditions_(std::move(levelset_boundary_conditions)) {
  /** Empty beside initializer list */
}

/**
 * @brief Performs all levelset halo updates in external boundaries.
 * @param node Node for which the halo update is done (indirect return).
 * @param buffer_type Identifier of the levelset buffer that should be updated.
 * @param loc Location of the boundary to update.
 */
void ExternalHaloManager::UpdateLevelsetExternal(
    Node &node, InterfaceBlockBufferType const buffer_type,
    BoundaryLocation const loc) const {
  levelset_boundary_conditions_[LTI(loc)]->UpdateLevelsetExternal(node,
                                                                  buffer_type);
}

/**
 * @brief Performs all interface tag halo updates in external boundaries.
 * @param interface_tags buffer for which the halo update is done (indirect
 * return).
 * @param loc Location of the boundary to update.
 */
void ExternalHaloManager::UpdateInterfaceTagExternal(
    std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
    BoundaryLocation const loc) const {
  levelset_boundary_conditions_[LTI(loc)]->UpdateInterfaceTagExternal(
      interface_tags);
}

/**
 * @brief Performs all interface tag halo updates in external boundaries.
 * @param Node for which the halo update is done (indirect return).
 * @param field_type Field identifier for material block buffers.
 * @param loc Location of the boundary to update.
 */
void ExternalHaloManager::UpdateMaterialExternal(
    Node &node, MaterialFieldType const field_type,
    BoundaryLocation const loc) const {
  material_boundary_conditions_[LTI(loc)]->UpdateMaterialExternal(node,
                                                                  field_type);
}
