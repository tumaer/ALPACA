//===------------------------ material_reader.h ---------------------------===//
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
#ifndef MATERIAL_READER_H
#define MATERIAL_READER_H

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "materials/equation_of_state_definitions.h"
#include "materials/material_property_definitions.h"
#include "materials/material_type_definitions.h"

/**
 * @brief Defines the class that provides access to the material data in the
 * input file. It serves as a proxy class for different material reader types
 * (xml,...) that only read the actual data. Here, consistency checks are done
 * that all read data are valid.
 */
class MaterialReader {

protected:
  // function to create the input tag of a given number of material indices
  std::string
  MaterialInputTag(std::vector<unsigned int> const &material_indices) const;

  // Functions that must be implemented by the derived classes
  virtual int DoReadNumberOfMaterials() const = 0;
  virtual std::string
  DoReadMaterialType(unsigned int const material_index) const = 0;
  virtual std::string
  DoReadEquationOfStateName(unsigned int const material_index) const = 0;
  virtual std::unordered_map<std::string, double>
  DoReadEquationOfStateData(unsigned int const material_index) const = 0;
  virtual double
  DoReadFixedValue(std::vector<unsigned int> const &material_indices,
                   MaterialProperty const property) const = 0;
  virtual std::string
  DoReadModelName(std::vector<unsigned int> const &material_indices,
                  MaterialProperty const property) const = 0;
  virtual std::unordered_map<std::string, double>
  DoReadModelData(std::vector<unsigned int> const &material_indices,
                  MaterialProperty const property) const = 0;

  // constructor can only be called from derived classes
  explicit MaterialReader() = default;

public:
  virtual ~MaterialReader() = default;
  MaterialReader(MaterialReader const &) = delete;
  MaterialReader &operator=(MaterialReader const &) = delete;
  MaterialReader(MaterialReader &&) = delete;
  MaterialReader &operator=(MaterialReader &&) = delete;

  // functions that return the different parameters
  TEST_VIRTUAL MaterialType ReadMaterialType(
      unsigned int const material_index, MaterialType const default_type) const;
  TEST_VIRTUAL unsigned int ReadNumberOfMaterials() const;
  TEST_VIRTUAL EquationOfStateName
  ReadEquationOfStateName(unsigned int const material_index) const;
  TEST_VIRTUAL std::unordered_map<std::string, double>
  ReadEquationOfStateData(unsigned int const material_index) const;
  TEST_VIRTUAL double
  ReadFixedValue(std::vector<unsigned int> const &material_indices,
                 MaterialProperty const property) const;
  MaterialPropertyModelName
  ReadModelName(std::vector<unsigned int> const &material_indices,
                MaterialProperty const property) const;
  std::unordered_map<std::string, double>
  ReadModelData(std::vector<unsigned int> const &material_indices,
                MaterialProperty const property) const;
};

#endif // MATERIAL_READER_H
