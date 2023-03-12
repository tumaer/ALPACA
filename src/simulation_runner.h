//===---------------------- simulation_runner.h ---------------------------===//
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
#ifndef SIMULATION_RUNNER_H
#define SIMULATION_RUNNER_H

#include "input_output/log_writer/log_writer.h"
#include "input_output/log_writer/logging.h"

#include "instantiation/halo_manager/instantiation_external_halo_manager.h"
#include "instantiation/halo_manager/instantiation_halo_manager.h"
#include "instantiation/halo_manager/instantiation_internal_halo_manager.h"
#include "instantiation/input_output/instantiation_input_output_manager.h"
#include "instantiation/input_output/instantiation_output_writer.h"
#include "instantiation/input_output/instantiation_restart_manager.h"
#include "instantiation/instantiation_communication_manager.h"
#include "instantiation/instantiation_initial_condition.h"
#include "instantiation/instantiation_modular_algorithm_assembler.h"
#include "instantiation/instantiation_multiresolution.h"
#include "instantiation/instantiation_unit_handler.h"
#include "instantiation/materials/instantiation_material_manager.h"
#include "instantiation/topology/instantiation_topology_manager.h"
#include "instantiation/topology/instantiation_tree.h"

#include "communication/mpi_utilities.h"
#include "user_specifications/space_filling_curve_settings.h"

namespace Simulation {

/**
 * @brief Run simulation of ALPACA.
 * @param input_reader Reader that is used to provide user-information from an
 * input file.
 */
void Run(InputReader const &input_reader) {
  LogWriter &logger = LogWriter::Instance();

  auto const input_file = input_reader.GetInputFile();

  logger.LogMessage("Using inputfile : " + input_file.string());

  auto const output_folder = InputOutput::CreateOutputBaseFolder(input_file);
  logger.SetLogfile(output_folder /
                    input_file.filename().replace_extension(".log"));
  logger.LogBreakLine();
  logger.Flush();

  Logging::LogCompiledSettings();

  // Create output_folder (base and set logfile)
  logger.LogMessage("Number of MPI ranks : " +
                    std::to_string(MpiUtilities::NumberOfRanks()));
  logger.LogMessage(
      "Load balancing      : " +
      SpaceFillingCurveSettings::SpaceFillingCurveSelectionString());
  logger.LogBreakLine();
  logger.Flush();

  // Instance for dimensionalization and non-dimensionalization of variables
  UnitHandler const unit_handler(
      Instantiation::InstantiateUnitHandler(input_reader));
  logger.LogBreakLine();
  logger.Flush();
  // Instance for handling of material and material pairing data
  MaterialManager const material_manager(
      Instantiation::InstantiateMaterialManager(input_reader, unit_handler));
  logger.LogBreakLine();
  logger.Flush();
  // Instance for handling global node data (cannot be const due to changes in
  // loop)
  TopologyManager topology_manager(Instantiation::InstantiateTopologyManager(
      input_reader, material_manager));
  Tree tree(Instantiation::InstantiateTree(input_reader, topology_manager,
                                           unit_handler));
  logger.LogBreakLine();
  logger.Flush();
  // Instance that provides multi-resolution information
  Multiresolution const multiresolution(
      Instantiation::InstantiateMultiresolution(input_reader,
                                                topology_manager));
  logger.LogBreakLine();
  logger.Flush();
  // Instance to provide communication
  CommunicationManager communication_manager(
      Instantiation::InstantiateCommunicationManager(topology_manager));
  // Instances for handling boundary conditions (internal and external). The
  // external and internal halo managers are not instantiated inside the halo
  // manager due to delete move constructors of both classes. A creation of the
  // external halo manager inside not suitable due to the inclusion of the input
  // reader. The internal cannot be created inside due to constness of some
  // functions.
  ExternalHaloManager const external_halo_manager(
      Instantiation::InstantiateExternalHaloManager(input_reader, unit_handler,
                                                    material_manager));
  InternalHaloManager internal_halo_manager(
      Instantiation::InstantiateInternalHaloManager(
          topology_manager, tree, communication_manager, material_manager));
  HaloManager halo_manager(Instantiation::InstantiateHaloManager(
      topology_manager, tree, external_halo_manager, internal_halo_manager,
      communication_manager));
  logger.LogBreakLine();
  logger.Flush();
  // Instance to restart simulation from snapshot and write output files (cannot
  // be const due to vector eraseing inside)
  OutputWriter const output_writer(Instantiation::InstantiateOutputWriter(
      topology_manager, tree, material_manager, unit_handler));
  RestartManager const restart_manager(Instantiation::InstantiateRestartManager(
      topology_manager, tree, unit_handler));
  InputOutputManager input_output_manager(
      Instantiation::InstantiateInputOutputManager(
          input_reader, output_writer, restart_manager, unit_handler,
          output_folder));
  logger.LogBreakLine();
  logger.Flush();
  // Instance for handling the initial conditions of the simulation
  std::unique_ptr<InitialCondition> initial_condition(
      Instantiation::InstantiateInitialCondition(input_reader, topology_manager,
                                                 tree, material_manager,
                                                 unit_handler));
  logger.LogBreakLine();
  logger.Flush();
  // Instance for the whole computation loop
  ModularAlgorithmAssembler mr_based_algorithm(
      Instantiation::InstantiateModularAlgorithmAssembler(
          input_reader, topology_manager, tree, communication_manager,
          halo_manager, multiresolution, material_manager, input_output_manager,
          unit_handler));
  logger.LogBreakLine();
  logger.Flush();

  // Initialize the simulation
  mr_based_algorithm.Initialization(*initial_condition);
  // Delete the initial condition by setting it to null
  initial_condition = nullptr;
  // Start loop computation
  mr_based_algorithm.ComputeLoop();

  logger.LogBreakLine();
  logger.Flush();
}

} // namespace Simulation

#endif // SIMULATION_RUNNER_H
