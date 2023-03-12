//===---------- instantiation_modular_algorithm_assembler.h ---------------===//
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
#ifndef INSTANTIATION_MODULAR_ALGORITHM_ASSEMBLER_H
#define INSTANTIATION_MODULAR_ALGORITHM_ASSEMBLER_H

#include "input_output/input_reader.h"
#include "modular_algorithm_assembler.h"

/**
 * @brief Defines all instantiation functions required for the modular algorithm
 * assembler.
 */
namespace Instantiation {

// Factory function to compute the gravity
std::array<double, 3> GetGravity(SourceTermReader const &source_term_reader,
                                 UnitHandler const &unit_handler);

// Instantiation function for the modular algorithm assembler
ModularAlgorithmAssembler InstantiateModularAlgorithmAssembler(
    InputReader const &input_reader, TopologyManager &topology_manager,
    Tree &tree, CommunicationManager &communication_manager,
    HaloManager &halo_manager, Multiresolution const &multiresolution,
    MaterialManager const &material_manager,
    InputOutputManager &input_output_manager, UnitHandler const &unit_handler);
} // namespace Instantiation

#endif // INSTANTIATION_MODULAR_ALGORITHM_ASSEMBLER_H
