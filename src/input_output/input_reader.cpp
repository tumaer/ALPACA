//===----------------------- input_reader.cpp -----------------------------===//
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
#include "input_output/input_reader.h"
#include <filesystem>

/**
 * @brief Standard constructor for the input reader class.
 * @param input_filename Name of the input file that is used.
 * @param input_type Identifier of the input type.
 * @param material_reader Base class to the material reader.
 * @param boundary_condition_reader Base class to the boundary condition reader.
 * @param initial_condition_reader Base class to the initial condition reader.
 * @param multi_resolution_reader Base class to the multi resolution reader.
 * @param dimensionalization_reader Base class to the dimensionalization reader.
 * @param output_reader Base class to the output reader.
 * @param restart_reader Base class to the restart reader.
 * @param source_term_reader Base class to the source term reader.
 * @param time_control_reader Base class to the time control reader.
 */
InputReader::InputReader(
    std::string const &input_filename, InputType const input_type,
    std::unique_ptr<MaterialReader const> material_reader,
    std::unique_ptr<BoundaryConditionReader const> boundary_condition_reader,
    std::unique_ptr<InitialConditionReader const> initial_condition_reader,
    std::unique_ptr<MultiResolutionReader const> multi_resolution_reader,
    std::unique_ptr<DimensionalizationReader const> dimensionalization_reader,
    std::unique_ptr<OutputReader const> output_reader,
    std::unique_ptr<RestartReader const> restart_reader,
    std::unique_ptr<SourceTermReader const> source_term_reader,
    std::unique_ptr<TimeControlReader const> time_control_reader)
    : // Start initializer list
      input_filename_(input_filename), input_type_(input_type),
      material_reader_(std::move(material_reader)),
      boundary_condition_reader_(std::move(boundary_condition_reader)),
      initial_condition_reader_(std::move(initial_condition_reader)),
      multi_resolution_reader_(std::move(multi_resolution_reader)),
      dimensionalization_reader_(std::move(dimensionalization_reader)),
      output_reader_(std::move(output_reader)),
      restart_reader_(std::move(restart_reader)),
      source_term_reader_(std::move(source_term_reader)),
      time_control_reader_(std::move(time_control_reader)) {
  /** Empty besides initializer list */
}

/**
 * @brief Gives the filename used as an input.
 * @return input filename.
 */
std::filesystem::path InputReader::GetInputFile() const {
  return input_filename_;
}

/**
 * @brief Gives the input type.
 * @return Input type identifier.
 */
InputType InputReader::GetInputType() const { return input_type_; }

/**
 * @brief Gives the instance of the material reader.
 * @return reference to the material reader.
 */
MaterialReader const &InputReader::GetMaterialReader() const {
  return *material_reader_;
}

/**
 * @brief Gives the instance of the boundary condition reader.
 * @return reference to the boundary condition reader.
 */
BoundaryConditionReader const &InputReader::GetBoundaryConditionReader() const {
  return *boundary_condition_reader_;
}

/**
 * @brief Gives the instance of the initial condition reader.
 * @return reference to the initial condition reader.
 */
InitialConditionReader const &InputReader::GetInitialConditionReader() const {
  return *initial_condition_reader_;
}

/**
 * @brief Gives the instance of the multiresolution reader.
 * @return reference to the multiresolution reader.
 */
MultiResolutionReader const &InputReader::GetMultiResolutionReader() const {
  return *multi_resolution_reader_;
}

/**
 * @brief Gives the instance of the non-dimensionalizaion reader.
 * @return reference to the non-dimensionalizaion reader.
 */
DimensionalizationReader const &
InputReader::GetDimensionalizationReader() const {
  return *dimensionalization_reader_;
}

/**
 * @brief Gives the instance of the output reader.
 * @return reference to the output reader.
 */
OutputReader const &InputReader::GetOutputReader() const {
  return *output_reader_;
}

/**
 * @brief Gives the instance of the restart reader.
 * @return reference to the restart reader.
 */
RestartReader const &InputReader::GetRestartReader() const {
  return *restart_reader_;
}

/**
 * @brief Gives the instance of the source term reader.
 * @return reference to the source term reader.
 */
SourceTermReader const &InputReader::GetSourceTermReader() const {
  return *source_term_reader_;
}

/**
 * @brief Gives the instance of the time control reader.
 * @return reference to the time control reader.
 */
TimeControlReader const &InputReader::GetTimeControlReader() const {
  return *time_control_reader_;
}
