//===---------------------- source_term_reader.h --------------------------===//
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
#ifndef SOURCE_TERM_READER_H
#define SOURCE_TERM_READER_H

#include <array>

#include "enums/direction_definition.h"

/**
 * @brief Defines the class that provides access to the source term data in the
 * input file. It serves as a proxy class for different source term reader types
 * (xml,...) that only read the actual data. Here, consistency checks are done
 * that all read data are valid.
 */
class SourceTermReader {

protected:
  // constructor can only be called from derived classes
  explicit SourceTermReader() = default;

  // Functions that must be implemented by the derived classes
  virtual double DoReadGravity(Direction const direction) const = 0;

public:
  virtual ~SourceTermReader() = default;
  SourceTermReader(SourceTermReader const &) = delete;
  SourceTermReader &operator=(SourceTermReader const &) = delete;
  SourceTermReader(SourceTermReader &&) = delete;
  SourceTermReader &operator=(SourceTermReader &&) = delete;

  // return function of the reader
  TEST_VIRTUAL double ReadGravity(Direction const direction) const;
};

#endif // SOURCE_TERM_READER_H
