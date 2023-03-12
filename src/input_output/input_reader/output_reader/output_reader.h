//===------------------------- output_reader.h ----------------------------===//
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
#ifndef OUTPUT_READER_H
#define OUTPUT_READER_H

#include <string>
#include <vector>

#include "input_output/output_writer/output_definitions.h"

/**
 * @brief Defines the class that provides access to the output data in the input
 * file. It serves as a proxy class for different output reader types (xml,...)
 * that only read the actual data. Here, consistency checks are done that all
 * read data are valid.
 */
class OutputReader {

protected:
  // constructor can only be called from derived classes
  explicit OutputReader() = default;

  // Functions that must be implemented by the derived classes
  virtual double DoReadTimeNamingFactor() const = 0;
  virtual std::string
  DoReadOutputTimesType(OutputType const output_type) const = 0;
  virtual double DoReadOutputInterval(OutputType const output_type) const = 0;
  virtual std::vector<double>
  DoReadOutputTimeStamps(OutputType const output_type) const = 0;

public:
  virtual ~OutputReader() = default;
  OutputReader(OutputReader const &) = delete;
  OutputReader &operator=(OutputReader const &) = delete;
  OutputReader(OutputReader &&) = delete;
  OutputReader &operator=(OutputReader &&) = delete;

  // Return functions proivded by the reader
  TEST_VIRTUAL double ReadTimeNamingFactor() const;
  TEST_VIRTUAL OutputTimesType
  ReadOutputTimesType(OutputType const output_type) const;
  TEST_VIRTUAL double ReadOutputInterval(OutputType const output_type) const;
  std::vector<double> ReadOutputTimeStamps(OutputType const output_type) const;
};

#endif // OUTPUT_READER_H
