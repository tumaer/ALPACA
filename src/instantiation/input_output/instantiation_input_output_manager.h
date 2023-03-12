//===-------------- instantiation_input_output_manager.h ------------------===//
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
#ifndef INSTANTIATION_INPUT_OUTPUT_MANAGER_H
#define INSTANTIATION_INPUT_OUTPUT_MANAGER_H

#include <filesystem>
#include <vector>

#include "input_output/input_output_manager.h"
#include "input_output/input_reader.h"

/**
 * @brief Defines all instantiation functions required for the input-output
 * manager.
 */
namespace Instantiation {

// factory functions for the input output manager
std::vector<double>
ComputeOutputTimes(OutputReader const &output_reader,
                   TimeControlReader const &time_control_reader,
                   UnitHandler const &unit_handler,
                   OutputType const output_type);
std::vector<double>
ComputeSnapshotTimes(RestartReader const &restart_reader,
                     TimeControlReader const &time_control_reader,
                     UnitHandler const &unit_handler);

// Instantiation function for the input_output manager
InputOutputManager InstantiateInputOutputManager(
    InputReader const &input_reader, OutputWriter const &output_writer,
    RestartManager const &restart_manager, UnitHandler const &unit_handler,
    std::filesystem::path base_output_folder);
} // namespace Instantiation

#endif // INSTANTIATION_INPUT_OUTPUT_MANAGER_H
