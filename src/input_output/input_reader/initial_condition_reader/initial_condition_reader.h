//===------------------ initial_condition_reader.h ------------------------===//
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
#ifndef INITIAL_CONDITION_READER_H
#define INITIAL_CONDITION_READER_H

#include <array>
#include <string>
#include <tuple>
#include <vector>

#include "initial_condition/levelset_initializer_definitions.h"
#include "initial_condition/parametric_variable.h"

/**
 * @brief Defines the class that provides access to the initial condition data
 * in the input file. It serves as a proxy class for different initial condition
 * reader types (xml,...) that only read the actual data. Here, consistency
 * checks are done that all read data are valid.
 */
class InitialConditionReader {

protected:
  // constructor can only be called from derived classes
  explicit InitialConditionReader() = default;

  // Functions that must be implemented by the derived classes
  virtual std::string
  DoReadMaterialInitialConditions(unsigned int const material_index) const = 0;
  virtual std::string
  DoReadLevelsetInitializerType(unsigned int const material_index) const = 0;
  virtual std::string
  DoReadLevelsetInitializerInput(unsigned int const levelset_index) const = 0;
  virtual std::vector<std::tuple<std::string, double, double, std::uint64_t>>
  DoReadParametricLevelsetInitializerVariables(
      unsigned int const levelset_index) const = 0;
  virtual std::array<double, 3>
  DoReadParametricLevelsetInitializerReferencePoint(
      unsigned int const levelset_index) const = 0;
  virtual std::vector<std::array<double, 6>>
  DoReadLevelsetInitializerBoundingBoxes(
      unsigned int const material_index) const = 0;

public:
  virtual ~InitialConditionReader() = default;
  InitialConditionReader(InitialConditionReader const &) = delete;
  InitialConditionReader &operator=(InitialConditionReader const &) = delete;
  InitialConditionReader(InitialConditionReader &&) = delete;
  InitialConditionReader &operator=(InitialConditionReader &&) = delete;

  // return functions of the reader class
  TEST_VIRTUAL std::string
  ReadMaterialInitialConditions(unsigned int const material_index) const;
  TEST_VIRTUAL LevelsetInitializerType
  ReadLevelsetInitializerType(unsigned int const levelset_index,
                              LevelsetInitializerType const default_type) const;
  TEST_VIRTUAL std::string
  ReadLevelsetInitializerInput(unsigned int const levelset_index) const;
  std::vector<ParametricVariable> ReadParametricLevelsetInitializerVariables(
      unsigned int const levelset_index) const;
  std::array<double, 3> ReadParametricLevelsetInitializerReferencePoint(
      unsigned int const levelset_index) const;
  TEST_VIRTUAL std::vector<std::array<double, 6>>
  ReadLevelsetInitializerBoundingBoxes(unsigned int const material_index) const;
};

#endif // INITIAL_CONDITION_READER_H
