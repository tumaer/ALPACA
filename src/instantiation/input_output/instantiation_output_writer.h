//===---------------- instantiation_output_writer.h -----------------------===//
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
#ifndef INSTANTIATION_OUTPUT_WRITER_H
#define INSTANTIATION_OUTPUT_WRITER_H

#include <memory>
#include <vector>

#include "input_output/output_writer.h"

/**
 * @brief Defines all instantiation functions required for the output writer.
 */
namespace Instantiation {

// factory functions for the output writer
std::unique_ptr<MeshGenerator const>
GetStandardMeshGenerator(UnitHandler const &unit_handler,
                         TopologyManager const &topology, Tree const &flower,
                         double const node_size_on_level_zero);
std::vector<std::unique_ptr<OutputQuantity const>>
GetMaterialOutputQuantities(UnitHandler const &unit_handler,
                            MaterialManager const &material_manager);
std::vector<std::unique_ptr<OutputQuantity const>>
GetInterfaceOutputQuantities(UnitHandler const &unit_handler,
                             MaterialManager const &material_manager);

// Instantiation function for the input_output manager
OutputWriter InstantiateOutputWriter(TopologyManager &topology_manager,
                                     Tree &tree,
                                     MaterialManager const &material_manager,
                                     UnitHandler const &unit_handler);
} // namespace Instantiation

#endif // INSTANTIATION_OUTPUT_WRITER_H
