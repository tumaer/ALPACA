//===-------------------- multi_resolution_reader.h -----------------------===//
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
#ifndef MULTI_RESOLUTION_READER_H
#define MULTI_RESOLUTION_READER_H

#include "enums/direction_definition.h"
#include <array>
#include <vector>

/**
 * @brief Defines the class that provides access to the multiresolution data in
 * the input file. It serves as a proxy class for different multiresolution
 * reader types (xml,...) that only read the actual data. Here, consistency
 * checks are done that all read data are valid.
 */
class MultiResolutionReader {

protected:
  // constructor can only be called from derived classes
  explicit MultiResolutionReader() = default;

  // Functions that must be implemented by the derived classes
  virtual double DoReadNodeSizeOnLevelZero() const = 0;
  virtual int DoReadNumberOfNodes(Direction const direction) const = 0;
  virtual int DoReadMaximumLevel() const = 0;
  virtual double DoReadEpsilonReference() const = 0;
  virtual int DoReadEpsilonLevelReference() const = 0;

public:
  virtual ~MultiResolutionReader() = default;
  MultiResolutionReader(MultiResolutionReader const &) = delete;
  MultiResolutionReader &operator=(MultiResolutionReader const &) = delete;
  MultiResolutionReader(MultiResolutionReader &&) = delete;
  MultiResolutionReader &operator=(MultiResolutionReader &&) = delete;

  // Function to return values with additional checks
  TEST_VIRTUAL double ReadNodeSizeOnLevelZero() const;
  TEST_VIRTUAL unsigned int ReadNumberOfNodes(Direction const direction) const;
  TEST_VIRTUAL unsigned int ReadMaximumLevel() const;
  TEST_VIRTUAL double ReadEpsilonReference() const;
  TEST_VIRTUAL unsigned int ReadEpsilonLevelReference() const;
};

#endif // MULTI_RESOLUTION_READER_H
