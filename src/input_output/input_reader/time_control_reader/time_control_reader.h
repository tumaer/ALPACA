//===--------------------- time_control_reader.h --------------------------===//
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
#ifndef TIME_CONTROL_READER_H
#define TIME_CONTROL_READER_H

/**
 * @brief Defines the class that provides access to the time control data in the
 * input file. It serves as a proxy class for different time control reader
 * types (xml,...) that only read the actual data. Here, consistency checks are
 * done that all read data are valid.
 */
class TimeControlReader {

protected:
  // Functions that must be implemented by the derived classes
  virtual double DoReadStartTime() const = 0;
  virtual double DoReadEndTime() const = 0;
  virtual double DoReadCFLNumber() const = 0;

  // constructor can only be called from derived classes
  explicit TimeControlReader() = default;

public:
  virtual ~TimeControlReader() = default;
  TimeControlReader(TimeControlReader const &) = default;
  TimeControlReader &operator=(TimeControlReader const &) = delete;
  TimeControlReader(TimeControlReader &&) = default;
  TimeControlReader &operator=(TimeControlReader &&) = delete;

  // Return functions with prepreparation of the data
  TEST_VIRTUAL double ReadStartTime() const;
  TEST_VIRTUAL double ReadEndTime() const;
  TEST_VIRTUAL double ReadCFLNumber() const;
};

#endif // TIME_CONTROL_READER_H
