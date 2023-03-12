//===--------------------- material_manager.h -----------------------------===//
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
#ifndef MATERIAL_MANAGER_H
#define MATERIAL_MANAGER_H

#include <memory>
#include <vector>

#include "materials/material.h"
#include "materials/material_definitions.h"
#include "materials/material_pairing.h"
#include "materials/material_type_definitions.h"
#include "unit_handler.h"

/**
 * @brief The MaterialManager class provides access to all materials and
 * material pairings present in the current simulation and forwards the
 * appropriate object to the caller. The MaterialManager does not change any
 * data. It provides the functionality to map material names and indices to the
 * correct material or pairing class.
 */
class MaterialManager {
  // Vector with all materials
  std::vector<std::tuple<MaterialType, Material>> const materials_;
  // Vector with all material pairings
  std::vector<MaterialPairing> const material_pairings_;

  // offset vector to provide proper mapping from material names to its material
  // pairings
  std::vector<int> const pairing_offset_;

  // local function to map a pairing of two materials to  corresponding vector
  // index of the material pairings
  unsigned int MapPairingToIndex(MaterialName const first_material,
                                 MaterialName const second_material) const;
  // factory function to generate the pairing offset vector (allows constness of
  // it)
  std::vector<int>
  GenerateMaterialPairingOffset(std::size_t const number_of_materials) const;

public:
  MaterialManager() = delete;
  explicit MaterialManager(
      std::vector<std::tuple<MaterialType, Material>> materials,
      std::vector<MaterialPairing> material_pairings);
  ~MaterialManager() = default;
  MaterialManager(MaterialManager const &) = delete;
  MaterialManager &operator=(MaterialManager const &) = delete;
  MaterialManager(MaterialManager &&) = delete;
  MaterialManager &operator=(MaterialManager &&) = delete;

  // Get the number of all materials contained in the simulation
  std::size_t GetNumberOfMaterials() const;

  // Get all materials used in the simulation
  std::vector<MaterialName> GetMaterialNames() const;

  // provides a material for a given index
  Material const &GetMaterial(std::size_t const index) const;

  // provides a material for a given identifier
  Material const &GetMaterial(MaterialName const material) const;

  // provides the material type for a given identifier
  MaterialType GetMaterialType(MaterialName const material) const;

  // provides the pairing of two material identifier
  MaterialPairing const &
  GetMaterialPairing(MaterialName const first_material,
                     MaterialName const second_material) const;

  // checks whether the given material is a solid boundary
  bool IsSolidBoundary(MaterialName const material) const;
};

#endif // MATERIAL_MANAGER_H
