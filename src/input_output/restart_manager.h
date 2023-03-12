//===------------------------- restart_manager.h --------------------------===//
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
#ifndef RESTART_MANAGER_H
#define RESTART_MANAGER_H

#include <string>

#include "input_output/hdf5/hdf5_manager.h"
#include "topology/topology_manager.h"
#include "topology/tree.h"
#include "unit_handler.h"

/**
 * @brief The RestartManager class handles the writing and reading of restart
 * snapshot files. The restart files are written in HDF5 Format (without XDMF
 * file). No mesh information is used for the restart file. Only data relevant
 * for setting up the simulations are required. The restart file cannot be
 * visualized in ParaView.
 */
class RestartManager {
  // Member variables to get topology information on current and all ranks
  Tree &tree_;
  TopologyManager &topology_;
  // Instance for logging to terminal/file
  LogWriter &logger_;
  // Instance to provide access to hdf5 functionality
  Hdf5Manager &hdf5_manager_;

  unsigned int const maximum_level_;
  double const length_reference_;
  double const density_reference_;
  double const temperature_reference_;
  double const velocity_reference_;

public:
  RestartManager() = delete;
  explicit RestartManager(UnitHandler const &unit_handler,
                          TopologyManager &topology_manager, Tree &tree,
                          unsigned int const maximum_level);
  ~RestartManager() = default;
  RestartManager(RestartManager const &) = delete;
  RestartManager &operator=(RestartManager const &) = delete;
  RestartManager(RestartManager &&) = delete;
  RestartManager &operator=(RestartManager &&) = delete;

  // Functions to restore simulation or write restart file
  double RestoreSimulation(std::string const &restore_filename) const;
  std::string
  WriteRestartFile(double const timestep,
                   std::string const &filename_without_extension) const;
};

#endif // RESTART_MANAGER_H
