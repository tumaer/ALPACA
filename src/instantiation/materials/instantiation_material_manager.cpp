//===--------------- instantiation_material_manager.cpp -------------------===//
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
#include "instantiation/materials/instantiation_material_manager.h"

#include "instantiation/materials/instantiation_material.h"
#include "instantiation/materials/instantiation_material_pairing.h"

namespace Instantiation {

/**
 * @brief Instantiates the complete set of materials with the given input
 * reader.
 * @param material_reader Reader that provides access to the material data of
 * the input file.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @return The fully instantiated set of materials present in the simulation.
 */
std::vector<std::tuple<MaterialType, Material>>
InstantiateMaterials(MaterialReader const &material_reader,
                     UnitHandler const &unit_handler) {

  // read the number of materials present in the simulation
  unsigned int const number_of_materials(
      material_reader.ReadNumberOfMaterials());

  // Check if too much materials have been specified
  if (number_of_materials + 1 >= MTI(MaterialName::MaterialOutOfBounds)) {
    throw std::invalid_argument(
        "This number of materials is not currently not supported!");
  }

  // declare vector that is returned and reserve enough memory
  std::vector<std::tuple<MaterialType, Material>> materials;
  materials.reserve(number_of_materials);
  // Instantiate each material individually
  for (unsigned int mat_index = 0; mat_index < number_of_materials;
       mat_index++) {
    materials.push_back(
        InstantiateMaterial(mat_index, material_reader, unit_handler));
  }

  // return the created vector
  return materials;
}

/**
 * @brief Instantiates the complete set of material pairings with the given
 * input reader.
 * @param material_reader Reader that provides access to the material data of
 * the input file.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @return The fully instantiated set of material pairings present in the
 * simulation.
 */
std::vector<MaterialPairing>
InstantiateMaterialPairings(MaterialReader const &material_reader,
                            UnitHandler const &unit_handler) {
  // First read the number of all materials
  unsigned int const number_of_materials =
      material_reader.ReadNumberOfMaterials();
  // Check if too much materials have been specified
  if (number_of_materials + 1 >= MTI(MaterialName::MaterialOutOfBounds)) {
    throw std::invalid_argument(
        "This number of materials is currently not supported!");
  }

  // Generate all combinations of indices
  std::vector<std::vector<unsigned int>> pairing_indices(
      GetMaterialPairingIndices(number_of_materials));

  // Check if too much material pairing combinations have been specified
  if (pairing_indices.size() + 1 >=
      MPTI(MaterialPairingName::MaterialPairingOutOfBounds)) {
    throw std::invalid_argument(
        "This number of material pairings is currently not supported!");
  }

  // declare vector that is returned
  std::vector<MaterialPairing> material_pairings;
  // Definition and Initialization of the materialPairing vector depending on
  // Single or Multiphase ( for later correct access ) NOTE: This if else can be
  // removed if the general multimaterial approach is incorporated in the
  // framework or the call to fill surface tension coefficient
  //       in capillary forces calculator is removed
  // Single material
  if (pairing_indices.size() == 0) {
    material_pairings.reserve(1);
    material_pairings.push_back(MaterialPairing());
  }
  // Multimaterial
  else {
    // Reserve enough entries in the material pairing vector
    material_pairings.reserve(pairing_indices.size());
    // Instantiate each pairing individually
    for (auto const &indices : pairing_indices) {
      material_pairings.push_back(
          InstantiateMaterialPairing(indices, material_reader, unit_handler));
    }
  }

  // return the created vector
  return material_pairings;
}

/**
 * @brief Instantiates the complete material manager class with the given input
 * reader.
 * @param input_reader Reader that provides access to the full data of the input
 * file.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @return The fully instantiated MaterialManager class.
 */
MaterialManager InstantiateMaterialManager(InputReader const &input_reader,
                                           UnitHandler const &unit_handler) {

  // First create and then move for proper logging inside (correct order)
  std::vector<std::tuple<MaterialType, Material>> materials(
      InstantiateMaterials(input_reader.GetMaterialReader(), unit_handler));
  std::vector<MaterialPairing> material_pairings(InstantiateMaterialPairings(
      input_reader.GetMaterialReader(), unit_handler));

  // Log a final empty line
  LogWriter &logger = LogWriter::Instance();
  logger.LogMessage(" ");

  return MaterialManager(std::move(materials), std::move(material_pairings));
}
} // namespace Instantiation
