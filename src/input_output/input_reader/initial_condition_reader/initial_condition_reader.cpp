//===----------------- initial_condition_reader.cpp -----------------------===//
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
#include "input_output/input_reader/initial_condition_reader/initial_condition_reader.h"

/**
 * @brief Reads the initial condition expression string for a given material
 * index.
 * @param material_index Index of the material for which the data is read.
 * @return Expression string for the material initial condition.
 */
std::string InitialConditionReader::ReadMaterialInitialConditions(
    unsigned int const material_index) const {
  return DoReadMaterialInitialConditions(material_index);
}

/**
 * @brief Gives the levelset initialization type name for a given levelset
 * index.
 * @param levelset_index The index of the levelset field (index start: 1).
 * @return Identifier of the levelset initialization type name.
 */
LevelsetInitializerType InitialConditionReader::ReadLevelsetInitializerType(
    unsigned int const levelset_index,
    LevelsetInitializerType const default_type) const {
  std::string const initializer_type(
      DoReadLevelsetInitializerType(levelset_index));
  return initializer_type.empty()
             ? default_type
             : StringToLevelsetInitializerType(initializer_type);
}

/**
 * @brief Gives the levelset initialization input for a given levelset index.
 * @param levelset_index The index of the levelset field (index start: 1).
 * @return The levelset initializer input string.
 */
std::string InitialConditionReader::ReadLevelsetInitializerInput(
    unsigned int const levelset_index) const {
  std::string const initializer_input(
      DoReadLevelsetInitializerInput(levelset_index));
  if (initializer_input.empty()) {
    throw std::invalid_argument(
        "The input for the levelset initial condition must be non-empty!");
  }
  return initializer_input;
}

/**
 * @brief Gives the variables for the parametric initializer.
 * @param levelset_index The index of the levelset (index start: 1).
 * @return The parametric variables.
 */
std::vector<ParametricVariable>
InitialConditionReader::ReadParametricLevelsetInitializerVariables(
    unsigned int const levelset_index) const {
  std::vector<std::tuple<std::string, double, double, std::uint64_t>> const
      parametric_variable_data(
          DoReadParametricLevelsetInitializerVariables(levelset_index));
  std::vector<ParametricVariable> parametric_variables(
      parametric_variable_data.size());
  std::transform(std::cbegin(parametric_variable_data),
                 std::cend(parametric_variable_data),
                 std::begin(parametric_variables), [](auto const &data_map) {
                   return ParametricVariable(
                       std::get<0>(data_map), std::get<1>(data_map),
                       std::get<2>(data_map), std::get<3>(data_map));
                 });
  return parametric_variables;
}

/**
 * @brief Gives the levelset bounding boxes for a given single levelset index.
 * @param levelset_index The index of the levelset (index start: 1).
 * @return Bounding boxes
 */
std::array<double, 3>
InitialConditionReader::ReadParametricLevelsetInitializerReferencePoint(
    unsigned int const levelset_index) const {
  return DoReadParametricLevelsetInitializerReferencePoint(levelset_index);
}

/**
 * @brief Gives the levelset bounding boxes for a given single levelset index.
 * @param levelset_index The index of the levelset (index start: 1).
 * @return Bounding boxes
 */
std::vector<std::array<double, 6>>
InitialConditionReader::ReadLevelsetInitializerBoundingBoxes(
    unsigned int const levelset_index) const {
  return DoReadLevelsetInitializerBoundingBoxes(levelset_index);
}
