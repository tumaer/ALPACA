//===------------------------- restart_reader.h ---------------------------===//
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
#ifndef RESTART_READER_H
#define RESTART_READER_H

#include "input_output/restart_manager/restart_definitions.h"
#include <string>
#include <vector>

/**
 * @brief Defines the class that provides access to the restart data in the
 * input file. It serves as a proxy class for different restart reader types
 * (xml,...) that only read the actual data. Here, consistency checks are done
 * that all read data are valid.
 */
class RestartReader {

protected:
  // constructor can only be called from derived classes
  explicit RestartReader() = default;

  // Functions that must be implemented by the derived classes
  virtual std::string DoReadRestoreMode() const = 0;
  virtual std::string DoReadRestoreFilename() const = 0;
  virtual std::string DoReadSnapshotTimesType() const = 0;
  virtual int DoReadSnapshotIntervalsToKeep() const = 0;
  virtual int DoReadSnapshotInterval() const = 0;
  virtual std::vector<double> DoReadSnapshotTimeStamps() const = 0;

public:
  virtual ~RestartReader() = default;
  RestartReader(RestartReader const &) = delete;
  RestartReader &operator=(RestartReader const &) = delete;
  RestartReader(RestartReader &&) = delete;
  RestartReader &operator=(RestartReader &&) = delete;

  // Return functions wiht included consistency checks
  TEST_VIRTUAL RestoreMode ReadRestoreMode() const;
  std::string ReadRestoreFilename() const;
  TEST_VIRTUAL SnapshotTimesType ReadSnapshotTimesType() const;
  unsigned int ReadSnapshotIntervalsToKeep() const;
  unsigned int ReadSnapshotInterval() const;
  TEST_VIRTUAL std::vector<double> ReadSnapshotTimeStamps() const;
};

#endif // RESTART_READER_H
