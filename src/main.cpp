//===---------------------------- main.cpp --------------------------------===//
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
#include "input_output/input_output_manager.h"
#include <fenv_wrapper.h> // Floating-Point raising exceptions.
#include <filesystem>
#include <mpi.h>

#include "communication/mpi_utilities.h"
#include "instantiation/input_output/instantiation_input_reader.h"
#include "instantiation/input_output/instantiation_log_writer.h"
#include "simulation_runner.h"

/**
 * @brief Starting function of ALPACA, called from the operating system.
 * @param argc Argument count set by your operating system
 * @param argv Input arguments, MPI settings, openMP settings, etc.
 * @return Zero if program finished correctly.
 * @note Please note, due to the usage of openMP #pragmas valgrind and such
 * tools might produce inaccurate measurements. This passes valgrind tests with
 * disabled openMP without errors
 */
int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  // Triggers signals on floating point errors, i.e. prohibits quiet NaNs and
  // alike
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // NH separate scope for MPI.
  {
    LogWriter &logger =
        Instantiation::InstantiateLogWriter(MpiUtilities::MasterRank());

    // determine the name of the executable and write it to the logger
    std::string const executable_name(argv[0]);
    logger.LogMessage("Using executable: " + executable_name);
    logger.Flush();

    // determine the name of the input file (default: inputfile.xml)
    std::filesystem::path const input_file(argc > 1 ? argv[1]
                                                    : "inputfile.xml");
    // Instance to provide interface to the input file/data
    InputReader const input_reader(
        Instantiation::InstantiateInputReader(input_file));

    Simulation::Run(input_reader);
    logger.Flush();
  }

  MPI_Finalize();

  return 0;
}
