//===------------------- dimensionalization_reader.h ----------------------===//
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
#ifndef DIMENSIONALIZATION_READER_H
#define DIMENSIONALIZATION_READER_H

/**
 * @brief Defines the class that provides access to the dimensionalization data
 * in the input file. It serves as a proxy class for different
 * dimensionalization reader types (xml,...) that only read the actual data.
 *        Here, consistency checks are done that all read data are valid.
 */
class DimensionalizationReader {

protected:
  // constructor can only be called from derived classes
  explicit DimensionalizationReader() = default;

  // Functions that must be implemented by the derived classes
  virtual double DoReadReferenceLength() const = 0;
  virtual double DoReadReferenceVelocity() const = 0;
  virtual double DoReadReferenceDensity() const = 0;
  virtual double DoReadReferenceTemperature() const = 0;

public:
  virtual ~DimensionalizationReader() = default;
  DimensionalizationReader(DimensionalizationReader const &) = delete;
  DimensionalizationReader &
  operator=(DimensionalizationReader const &) = delete;
  DimensionalizationReader(DimensionalizationReader &&) = delete;
  DimensionalizationReader &operator=(DimensionalizationReader &&) = delete;

  // Return functions
  double ReadReferenceLength() const;
  double ReadReferenceVelocity() const;
  double ReadReferenceDensity() const;
  double ReadReferenceTemperature() const;
};

#endif // DIMENSIONALIZATION_READER_H
