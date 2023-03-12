//===---------- instantiation_modular_algorithm_assembler.cpp -------------===//
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
#include "instantiation/instantiation_modular_algorithm_assembler.h"

#include "user_specifications/compile_time_constants.h"
#include "utilities/string_operations.h"

namespace Instantiation {

/**
 * @brief Computes the gravity source term.
 * @param source_term_reader The reader providing access to all source term
 * information of the input file. param unit_handler Instance to provide
 * (non-)dimensionalization of values.
 * @return gravity in x-, y- and z-direction.
 */
std::array<double, 3> GetGravity(SourceTermReader const &source_term_reader,
                                 UnitHandler const &unit_handler) {
  // Initialize the gravity with zero
  std::array<double, 3> gravity = {0.0, 0.0, 0.0};

  // Only proceed if gravity is activated
  if constexpr (CC::GravityIsActive()) {
    gravity[0] = unit_handler.NonDimensionalizeValue(
        source_term_reader.ReadGravity(Direction::X), {UnitType::Length},
        {UnitType::Time, UnitType::Time});
    // For two and three dimensions
    if constexpr (CC::DIM() != Dimension::One) {
      gravity[1] = unit_handler.NonDimensionalizeValue(
          source_term_reader.ReadGravity(Direction::Y), {UnitType::Length},
          {UnitType::Time, UnitType::Time});
    }
    // Only for three dimensions
    if constexpr (CC::DIM() == Dimension::Three) {
      gravity[2] = unit_handler.NonDimensionalizeValue(
          source_term_reader.ReadGravity(Direction::Z), {UnitType::Length},
          {UnitType::Time, UnitType::Time});
    }
  }

  return gravity;
}

/**
 * @brief Get a vector containing all levels of the current simulation.
 * @param maximum_level Maximum level used for the simulation.
 * @return vector with all levels.
 */
std::vector<unsigned int> GetAllLevels(unsigned int const maximum_level) {
  std::vector<unsigned int> all_levels(
      maximum_level + 1); // Level zero needs to be counted as well
  std::iota(all_levels.begin(), all_levels.end(), 0);
  return all_levels;
}

/**
 * @brief Instantiates the complete modular algorithm assembler class with the
 * given input classes.
 * @param input_reader Reader that provides access to the full data of the input
 * file.
 * @param topology_manager Class providing global (on all ranks) node
 * information.
 * @param tree Tree class providing local (on current rank) node information.
 * @param communication_manager Calls providing communication handling between
 * different ranks.
 * @param multiresolution Instance to provide mutliresolution computations for
 * remeshing and so fourth.
 * @param material_manager Instance providing instantiated material data.
 * @param input_output_manager Instance to provide restart handling and output
 * writing.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @return The fully instantiated Modular algorithm assembler class.
 */
ModularAlgorithmAssembler InstantiateModularAlgorithmAssembler(
    InputReader const &input_reader, TopologyManager &topology_manager,
    Tree &tree, CommunicationManager &communication_manager,
    HaloManager &halo_manager, Multiresolution const &multiresolution,
    MaterialManager const &material_manager,
    InputOutputManager &input_output_manager, UnitHandler const &unit_handler) {

  // Get data that is logged
  double const start_time = unit_handler.NonDimensionalizeValue(
      input_reader.GetTimeControlReader().ReadStartTime(), UnitType::Time);
  double const end_time = unit_handler.NonDimensionalizeValue(
      input_reader.GetTimeControlReader().ReadEndTime(), UnitType::Time);
  double const cfl_number = input_reader.GetTimeControlReader().ReadCFLNumber();

  // Log data
  LogWriter &logger = LogWriter::Instance();
  logger.LogMessage(" ");
  logger.LogMessage("Time control data: ");
  logger.LogMessage(
      StringOperations::Indent(2) + "Start time: " +
      StringOperations::ToScientificNotationString(
          unit_handler.DimensionalizeValue(start_time, UnitType::Time), 9));
  logger.LogMessage(
      StringOperations::Indent(2) + "End time  : " +
      StringOperations::ToScientificNotationString(
          unit_handler.DimensionalizeValue(end_time, UnitType::Time), 9));
  logger.LogMessage(
      StringOperations::Indent(2) + "CFL number: " +
      StringOperations::ToScientificNotationString(cfl_number, 9));
  logger.LogMessage(" ");
  // Compute the cell size on maximum level
  unsigned int const maximum_level = topology_manager.GetMaximumLevel();
  // The MAA required the dimensionless cell size -> No dimensionalization
  // required
  double const cell_size_on_maximum_level = tree.GetNodeSizeOnLevelZero() /
                                            double(CC::ICX()) /
                                            double(1 << maximum_level);

  // initialize the algorithm assembler
  return ModularAlgorithmAssembler(
      start_time, end_time, cfl_number,
      GetGravity(input_reader.GetSourceTermReader(), unit_handler),
      GetAllLevels(maximum_level), cell_size_on_maximum_level, unit_handler,
      tree, topology_manager, halo_manager, communication_manager,
      multiresolution, material_manager, input_output_manager);
}
} // namespace Instantiation
