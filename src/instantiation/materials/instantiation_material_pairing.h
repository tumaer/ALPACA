//===----------------- instantiation_material_pairing.h -------------------===//
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
#ifndef INITIALIZATION_MATERIAL_PAIRING_H
#define INITIALIZATION_MATERIAL_PAIRING_H

#include "input_output/input_reader/material_reader/material_reader.h"
#include "materials/material_pairing.h"
#include <memory>

/**
 * @brief Defines all instantiation functions required for one single material
 * pairing.
 */
namespace Instantiation {

// initialize function for the surface tension coefficient model
std::unique_ptr<InterfaceParameterModel const>
InstantiateSurfaceTensionCoefficientModel(
    MaterialPropertyModelName const model_name,
    std::unordered_map<std::string, double> const &model_data,
    UnitHandler const &unit_handler);

// initialize function for the complete material pairing
MaterialPairing
InstantiateMaterialPairing(std::vector<unsigned int> const &material_indices,
                           MaterialReader const &material_reader,
                           UnitHandler const &unit_handler);

} // namespace Instantiation

#endif // INITIALIZATION_MATERIAL_PAIRING_H
