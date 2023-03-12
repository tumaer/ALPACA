//===--------------------- material_definitions.h -------------------------===//
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
#ifndef MATERIAL_DEFINITIONS_H
#define MATERIAL_DEFINITIONS_H

#include <type_traits>
#include <vector>

/**
 * @brief The MaterialName enum gives all materials a unique identifier to allow
 * automatic selection of correct materials as well as to distinguish treatment
 * in case of contact between materials, e.g. merging or interaction.
 * @note !IMPORTANT! Do NOT change the order of assignment from MaterialOne to
 * index 0, MaterialTwo to index 1, and so fourth. This is used to index the
 * vector in the material manager class and for assignment of proper material
 * pairings between two material. Adding of additional material by simply
 * increasing numbers of MaterialOutOfBounds and filling the gap in between.
 */
enum class MaterialName : unsigned short {
  MaterialOne = 0,
  MaterialTwo = 1,
  MaterialThree = 2,
  MaterialFour = 3,
  MaterialOutOfBounds = 4
};

/**
 * @brief Converts a Material identifier to a (C++11 standard compliant, i. e.
 * positive) sortabel integer. "MTI = Material To Index".
 * @param m The material identifier.
 * @return Integer to be used for sorting.
 */
constexpr std::underlying_type<MaterialName>::type MTI(MaterialName const m) {
  return static_cast<typename std::underlying_type<MaterialName>::type>(m);
}

/**
 * @brief Converts a Material index to its corresponding enum name "ITM = Index
 * to Material".
 * @param i The material index.
 * @return Material name for the index.
 */
constexpr MaterialName
ITM(typename std::underlying_type<MaterialName>::type i) {
  return static_cast<MaterialName>(i);
}

/**
 * @brief The MaterialPairingName enum gives all identifiers for a pairing
 * between two materials.
 * @note !IMPORTANT! Do NOT change the order of assignment from MaterialOneTwo
 * to index 0, MaterialOneThree to index 1, MaterialTwoThree to index 2 and so
 * fourth. This is used to index the vector in the material manager class and
 * provides a unique mapping function. If an additionall material is added
 * simply add all pairings in the order MaterialPairing[1:N-1][N] and increase
 * the out of bounds number.
 */
enum class MaterialPairingName : unsigned short {
  MaterialPairingOneTwo = 0,
  MaterialPairingOneThree = 1,
  MaterialPairingTwoThree = 2,
  MaterialPairingOneFour = 3,
  MaterialPairingTwoFour = 4,
  MaterialPairingThreeFour = 5,
  MaterialPairingOutOfBounds = 6
};

/**
 * @brief Converts a MaterialPairing identifier to a (C++11 standard compliant,
 * i. e. positive) sortabel integer. "MPTI = Material Pairing To Index".
 * @param mp The material pairing identifier.
 * @return Integer to be used for sorting.
 */
constexpr std::underlying_type<MaterialPairingName>::type
MPTI(MaterialPairingName const mp) {
  return static_cast<typename std::underlying_type<MaterialPairingName>::type>(
      mp);
}

/**
 * @brief Converts a MaterialPairing index to its corresponding enum name "ITMP
 * = Index to Material Pairing".
 * @param i The material pairing index.
 * @return Name of the material pairing.
 */
constexpr MaterialPairingName
ITMP(typename std::underlying_type<MaterialPairingName>::type i) {
  return static_cast<MaterialPairingName>(i);
}

/**
 * @brief Creates for the given number of materials all pairing indices.
 * @param number_of_materials Number of materials considered.
 * @return vector with pairing indices.
 *
 * @note Do not change the order here, since this reflects the same order as in
 * the enum class of the pairings.
 */
inline std::vector<std::vector<unsigned int>>
GetMaterialPairingIndices(unsigned int const number_of_materials) {
  // Generate all combinations of indices
  std::vector<std::vector<unsigned int>> pairing_indices;
  // Loop structure for creation
  for (unsigned int second_index = 1; second_index < number_of_materials;
       ++second_index) {
    for (unsigned int first_index = 0; first_index <= second_index - 1;
         ++first_index) {
      pairing_indices.push_back({first_index + 1, second_index + 1});
    }
  }

  return pairing_indices;
}

#endif // MATERIAL_DEFINITIONS_H
