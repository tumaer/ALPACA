//===------------------- boundary_condition_reader.h ----------------------===//
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
#ifndef BOUNDARY_CONDITION_READER_H
#define BOUNDARY_CONDITION_READER_H

#include <array>

#include "block_definitions/field_material_definitions.h"
#include "boundary_condition/boundary_specifications.h"
#include "enums/direction_definition.h"

/**
 * @brief Defines the class that provides access to the boundary condition data
 * in the input file. It serves as a proxy class for different boundary
 * condition reader types (xml,...) that only read the actual data. Here,
 * consistency checks are done that all read data are valid.
 */
class BoundaryConditionReader {

protected:
  // constructor can only be called from derived classes
  explicit BoundaryConditionReader() = default;

  // Functions that must be implemented by the derived classes
  virtual std::string
  DoReadMaterialBoundaryType(BoundaryLocation const location) const = 0;
  virtual std::string
  DoReadLevelSetBoundaryType(BoundaryLocation const location) const = 0;
  virtual double DoReadMaterialFixedValueBoundaryCondition(
      BoundaryLocation const location, std::string const &variable) const = 0;

public:
  virtual ~BoundaryConditionReader() = default;
  BoundaryConditionReader(BoundaryConditionReader const &) = delete;
  BoundaryConditionReader &operator=(BoundaryConditionReader const &) = delete;
  BoundaryConditionReader(BoundaryConditionReader &&) = delete;
  BoundaryConditionReader &operator=(BoundaryConditionReader &&) = delete;

  // Functions that must be implemented by the derived classes
  TEST_VIRTUAL MaterialBoundaryType
  ReadMaterialBoundaryType(BoundaryLocation const location) const;
  TEST_VIRTUAL LevelSetBoundaryType
  ReadLevelsetBoundaryType(BoundaryLocation const location) const;
  TEST_VIRTUAL std::array<double, MF::ANOP()>
  ReadMaterialFixedValueBoundaryConditions(
      BoundaryLocation const location) const;
};

#endif // BOUNDARY_CONDITION_READER_H
