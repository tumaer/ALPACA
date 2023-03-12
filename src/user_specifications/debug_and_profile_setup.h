//===-------------------- debug_and_profile_setup.h -----------------------===//
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
#ifndef DEBUG_PROFILE_SETUP_H
#define DEBUG_PROFILE_SETUP_H

/**
 * @brief Small class to control compile-time settings related with additional
 * logging and output for profiling or debugging reasons.
 */
class DebugProfileSetup {

  static constexpr bool debug_logging_ =
      false; // Writes debug information to the logger
  static constexpr bool debug_output_ =
      false; // Writes the debug output (mesh and simulation data) to file
  static constexpr bool profiling_ = false;

public:
  /**
   * @brief Gives a bool indicating if any debugging option is set.
   * @return True if any debugging option is enabled. False otherwise.
   */
  static constexpr bool Debug() { return debug_logging_ || debug_output_; }

  /**
   * @brief Gives a bool to decide if additional debugging information should be
   * provided. This function is not about debugging information for e.g. GNU
   * Debugger, but to print additional information during a debug run or to
   * create additional information after separate functions. All information is
   *        written to the logger class.
   *        $Convenience function. The compiler should recognize unset debug if
   * clauses and take them out in release-builds (hopefully).$
   * @return Debug logging decision for the current build.
   */
  static constexpr bool DebugLog() { return debug_logging_; }

  /**
   * @brief Gives a bool to decide if additional debug output files are to be
   * written out. The output is the actual mesh and simulation data of all cells
   * that are written into the output format specified through the input file
   * (e.g., XDMF/HDF5)
   * @return Debug output decision for the current build.
   */
  static constexpr bool DebugOutput() { return debug_output_; }

  /**
   * @brief Give a bool to indicate whether or not profiling (timing) runs are
   * beeing run.
   * @return Profiling decision.
   */
  static constexpr bool Profile() { return profiling_; }
};

using DP = DebugProfileSetup;

#endif // DEBUG_PROFILE_SETUP_H
