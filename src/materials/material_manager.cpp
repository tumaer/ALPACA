//===----------------------- material_manager.cpp -------------------------===//
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
#include "materials/material_manager.h"

#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "utilities/string_operations.h"
#include <stdexcept>

/**
 * @brief Sets up a MaterialManager for the given pure materials and material
 * pairings.
 * @param materials All already initialized materials.
 * @param material_pairing_data All already initialized material pairings.
 */
MaterialManager::MaterialManager(
    std::vector<std::tuple<MaterialType, Material>> materials,
    std::vector<MaterialPairing> material_pairings)
    : // Start initializer list
      materials_(std::move(materials)),
      material_pairings_(std::move(material_pairings)),
      pairing_offset_(GenerateMaterialPairingOffset(materials_.size())) {

  // Instantiation of material sign capsule for two-material computations
  // AB 20-03-30: This concept could be used everywhere if material manager is
  // present in class and must be implemented in case of
  //              multi-material implementation
  MaterialSignCapsule(ITM(0), ITM(materials_.size() - 1));
}

/**
 * @brief Maps the pairing of two materials to the correct position of the
 * vector using the material enum class index and a pairing offset.
 * @param first_material, second_material materials for which the pairing should
 * be obtained.
 * @note This mapping ensures that the same result is ensured regardless of the
 * order of both materials in the function call.
 */
unsigned int
MaterialManager::MapPairingToIndex(MaterialName const first_material,
                                   MaterialName const second_material) const {

  // NOTE: This part is required to have a fall back behavior in case capillary
  // forces are activated and only one material is defined.
  if (materials_.size() == 1) {
    return 0;
  }

#ifndef PERFORMANCE
  if (first_material == second_material) {
    throw std::runtime_error("Error! Material pairing can only be obtained for "
                             "two distinct materials!");
  }
#endif

  // Obtain the indices of the two materials in the materials vector (
  // implicitly done due to the underlying index of MaterialName enum )
  auto material_A = MTI(first_material);
  auto material_B = MTI(second_material);

  // return the mapping index ( conversion to unsigned int is safe, since the
  // max function always returns at least index 1 )
  return material_A + material_B +
         pairing_offset_[std::max(material_A, material_B) - 1];
}

/**
 * @brief Returns the total number of materials contained in the simulation.
 * @return Number of materials.
 */
std::size_t MaterialManager::GetNumberOfMaterials() const {
  return materials_.size();
}

/**
 * @brief Returns all materials contained in the simulation.
 * @return Vector with all material names.
 */
std::vector<MaterialName> MaterialManager::GetMaterialNames() const {
  // Declare vector that is returned
  std::vector<MaterialName> material_names;
  // fill vector with all names
  material_names.reserve(materials_.size());
  for (unsigned int material_index = 0; material_index < materials_.size();
       material_index++) {
    material_names.push_back(ITM(material_index));
  }

  return material_names;
}

/**
 * @brief Gives an instance to the material with the given index.
 * @param index Index of the material that should be returned.
 * @return Instance of the material.
 * @note No sanity check is done here that the index really exists. Safety is
 * ensured if the indices are obtained with the GetNumberOfMaterials() function.
 */
Material const &MaterialManager::GetMaterial(std::size_t const index) const {

  return std::get<1>(materials_[index]);
}

/**
 * @brief Gives an instance of the material with the given identifier.
 * @param material The material identifier.
 * @return Instance of the material.
 * @note No sanity check is done here, since in general all present material
 * names that are created should be checked by the MaterialName enum class
 * itself.
 */
Material const &
MaterialManager::GetMaterial(MaterialName const material) const {

  return std::get<1>(materials_[MTI(material)]);
}

/**
 * @brief Gives the material type of the material with the given identifier.
 * @param material The material identifier.
 * @return The material type.
 * @note No sanity check is done here, since in general all present material
 * names that are created should be checked by the MaterialName enum class
 * itself.
 */
MaterialType
MaterialManager::GetMaterialType(MaterialName const material) const {

  return std::get<0>(materials_[MTI(material)]);
}

/**
 * @brief Gives an instance of a material pairing for two different materials.
 * @param first_material, second_material Unique identifiers of the materials of
 * interest.
 * @return Instance of the material pairing.
 */
MaterialPairing const &
MaterialManager::GetMaterialPairing(MaterialName const first_material,
                                    MaterialName const second_material) const {
  return material_pairings_[MapPairingToIndex(first_material, second_material)];
}

/**
 * @brief
 * @return
 */
bool MaterialManager::IsSolidBoundary(MaterialName const material) const {
  return GetMaterialType(material) == MaterialType::SolidBoundary;
}

/**
 * @brief Generates the offset values for each material to provide proper
 * mapping of the material indices towards material pairing indices.
 * @param number_of_materials Number of materials to be considered.
 * @return Vector with the appropriate offset values for the given number of
 * materials.
 */
std::vector<int> MaterialManager::GenerateMaterialPairingOffset(
    size_t const number_of_materials) const {
  // declare vector
  std::vector<int> pairing_offset;
  // Definition and Initialization of the materialPairing vector depending on
  // Single or Multiphase ( for later correct access ) NOTE: This if else can be
  // removed if the general multimaterial approach is incorporated in the
  // framework. Single material
  if (number_of_materials == 1) {
    pairing_offset.reserve(1);
    pairing_offset.push_back(0);
  }
  // Multimaterial
  else {
    // Define the pairing offset for later mapping of the pairing to the
    // appropriate positions in vector
    pairing_offset.reserve(number_of_materials - 1);

    // The offset represents the increases of the index to access
    // material_pairings_ vector resulting from an increase in material numbers
    // compared to the previous number of materials. The pairing_offset is then
    // used together with the index of the materials enum class to access the
    // correct position in the material_pairings vector.
    // 2-materials: -1, 3-materials: -1, 4-materials: 0, 5-materials: 2;
    // 6-fludis: 5, 7-materials: 9, ...
    int offset = -1;
    pairing_offset.push_back(offset);
    for (size_t n = 3; n <= number_of_materials; n++) {
      // increment offset
      offset += n - 3;
      // add offset at correct position for number of materials
      pairing_offset.push_back(offset);
    }
  }

  // return the vector
  return pairing_offset;
}

// MaterialSignCapsule does not have a cpp file for definition
MaterialName MaterialSignCapsule::positive_material_;
MaterialName MaterialSignCapsule::negative_material_;
