//===---------------- boundary_condition_reader.cpp -----------------------===//
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
#include "input_output/input_reader/boundary_condition_reader/boundary_condition_reader.h"

/**
 * @brief Gives the boundary type of the materials at the given location.
 * @param location Boundary location for which the type is read (east, west,
 * north, south, top, bottom).
 * @return Material boundary type.
 */
MaterialBoundaryType BoundaryConditionReader::ReadMaterialBoundaryType(
    BoundaryLocation const location) const {
  return StringToMaterialBoundaryType(DoReadMaterialBoundaryType(location));
}

/**
 * @brief Gives the boundary type of the levelset at the given location.
 * @param location Boundary location for which the type is read (east, west,
 * north, south, top, bottom).
 * @return Levelset boundary type.
 */
LevelSetBoundaryType BoundaryConditionReader::ReadLevelsetBoundaryType(
    BoundaryLocation const location) const {
  return StringToLevelSetBoundaryType(DoReadLevelSetBoundaryType(location));
}

/**
 * @brief Reads the (reduced) set of fixed value prime states at a given
 * location.
 * @param location Boundary location for which the type is read (east, west,
 * north, south, top, bottom).
 * @return Reduced set of active prime states.
 */
std::array<double, MF::ANOP()>
BoundaryConditionReader::ReadMaterialFixedValueBoundaryConditions(
    BoundaryLocation const location) const {

  // declare the array for the fixed values
  std::array<double, MF::ANOP()> fixed_values;

  // read all prime state values from the input file
  for (PrimeState const prime : MF::ASOP()) {
    // Obtain the name of the prime state
    std::string const prime_name(std::string(MF::InputName(prime)));
    // First check if the variable should be read or not (empty -> no reading)
    if (prime_name.empty()) {
      fixed_values[PTI(prime)] = 0.0;
    } else {
      fixed_values[PTI(prime)] =
          DoReadMaterialFixedValueBoundaryCondition(location, prime_name);
    }
  }

  return fixed_values;
}
