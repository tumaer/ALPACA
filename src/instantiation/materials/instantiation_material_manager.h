//===----------------- instantiation_material_manager.h -------------------===//
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
#ifndef INITIALIZATION_MATERIAL_MANAGER_H
#define INITIALIZATION_MATERIAL_MANAGER_H

#include "input_output/input_reader.h"
#include "materials/material_manager.h"

/**
 * @brief Defines all instantiation functions required for the material manager.
 */
namespace Instantiation {

// Instantiation function fot the full set of materials
std::vector<std::tuple<MaterialType, Material>>
InstantiateMaterials(MaterialReader const &material_reader,
                     UnitHandler const &unit_handler);

// Instantiation function for the full set of material pairings
std::vector<MaterialPairing>
InstantiateMaterialPairings(MaterialReader const &material_reader,
                            UnitHandler const &unit_handler);

// Instantiation function for the complete material manager
MaterialManager InstantiateMaterialManager(InputReader const &input_reader,
                                           UnitHandler const &unit_handler);

} // namespace Instantiation

#endif // INITIALIZATION_MATERIAL_MANAGER_H
