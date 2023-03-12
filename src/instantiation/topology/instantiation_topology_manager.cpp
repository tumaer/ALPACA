//===------------- instantiation_topology_manager.cpp ---------------------===//
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
#include "instantiation/topology/instantiation_topology_manager.h"

#include "topology/id_periodic_information.h"

namespace Instantiation {

/**
 * @brief Computes the correct number of blocks on level zero used in the
 * simulation.
 * @param multi_resolution_reader Instance that provides access to the
 * multiresolution data in the input file.
 * @return Number of blocks in each direction.
 */
std::array<unsigned int, 3> GetNumberOfNodesOnLevelZero(
    MultiResolutionReader const &multi_resolution_reader) {
  // Initialize number of blocks with default values
  std::array<unsigned int, 3> number_of_nodes = {1, 1, 1};
  // Read each direction if required
  number_of_nodes[0] = multi_resolution_reader.ReadNumberOfNodes(Direction::X);
  // For two and three dimensions
  if constexpr (CC::DIM() != Dimension::One) {
    number_of_nodes[1] =
        multi_resolution_reader.ReadNumberOfNodes(Direction::Y);
  }
  // Only for three dimensions
  if constexpr (CC::DIM() == Dimension::Three) {
    number_of_nodes[2] =
        multi_resolution_reader.ReadNumberOfNodes(Direction::Z);
  }

  return number_of_nodes;
}

/**
 * @brief Checks whether the material boundary conditions are periodic or not
 * for a given direction.
 * @param boundary_condition_reader Instance that provide access to the boundary
 * condition data in the input file.
 * @param direction Direction for which periodicity is checked.
 * @return true if it is periodic, false else.
 * @note Function throws error if boundary types on opposite sides do not match.
 */
bool IsMaterialPeriodic(
    BoundaryConditionReader const &boundary_condition_reader,
    Direction const direction) {
  // Declare the boundary locations
  BoundaryLocation const first_location =
      direction == Direction::X   ? BoundaryLocation::East
      : direction == Direction::Y ? BoundaryLocation::North
                                  : BoundaryLocation::Bottom;
  BoundaryLocation const second_location =
      direction == Direction::X   ? BoundaryLocation::West
      : direction == Direction::Y ? BoundaryLocation::South
                                  : BoundaryLocation::Top;
  // Declare the two boundary conditions
  MaterialBoundaryType const first_boundary =
      boundary_condition_reader.ReadMaterialBoundaryType(first_location);
  MaterialBoundaryType const second_boundary =
      boundary_condition_reader.ReadMaterialBoundaryType(second_location);

  // check periodicity
  if (first_boundary == MaterialBoundaryType::Periodic ||
      second_boundary == MaterialBoundaryType::Periodic) {
    if (first_boundary != second_boundary) {
      throw std::invalid_argument(
          "Incorrect use of " + BoundaryLocationToString(first_location, true) +
          "-" + BoundaryLocationToString(second_location, true) +
          " periodic condition, both boundaries from the material must be "
          "periodic!");
    } else {
      return true;
    }
  } else {
    return false;
  }
}

/**
 * @brief Checks whether the levelset boundary conditions are periodic or not
 * for a given direction.
 * @param boundary_condition_reader Instance that provide access to the boundary
 * condition data in the input file.
 * @param direction Direction for which periodicity is checked.
 * @return true if it is periodic, false else.
 * @note Function throws error if boundary types on opposite sides do not match.
 */
bool IsLevelsetPeriodic(
    BoundaryConditionReader const &boundary_condition_reader,
    Direction const direction) {
  // Declare the boundary locations
  BoundaryLocation const first_location =
      direction == Direction::X   ? BoundaryLocation::East
      : direction == Direction::Y ? BoundaryLocation::North
                                  : BoundaryLocation::Bottom;
  BoundaryLocation const second_location =
      direction == Direction::X   ? BoundaryLocation::West
      : direction == Direction::Y ? BoundaryLocation::South
                                  : BoundaryLocation::Top;
  // Declare the two boundary conditions
  LevelSetBoundaryType const first_boundary =
      boundary_condition_reader.ReadLevelsetBoundaryType(first_location);
  LevelSetBoundaryType const second_boundary =
      boundary_condition_reader.ReadLevelsetBoundaryType(second_location);

  // check periodicity
  if (first_boundary == LevelSetBoundaryType::Periodic ||
      second_boundary == LevelSetBoundaryType::Periodic) {
    if (first_boundary != second_boundary) {
      throw std::invalid_argument(
          "Incorrect use of " + BoundaryLocationToString(first_location, true) +
          "-" + BoundaryLocationToString(second_location, true) +
          " periodic condition, both boundaries from the levelset must be "
          "periodic!");
    } else {
      return true;
    }
  } else {
    return false;
  }
}

/**
 * @brief Gives the value for the active periodic directions used in the
 * simulation.
 * @param boundary_condition_reader Reader to access the boundary condition data
 * of the input file.
 * @param material_manager material_manager Instance providing initialized
 * material data.
 * @return Value of active periodic directions.
 */
unsigned int GetActivePeriodicDirections(
    BoundaryConditionReader const &boundary_condition_reader,
    MaterialManager const &material_manager) {
  // Define the return value
  unsigned int active_periodic_locations = 0;

  // read the number of materials
  std::size_t const number_of_materials =
      material_manager.GetNumberOfMaterials();

  // Check the x-direction (levelset only if multi material)
  bool levelset_periodic =
      number_of_materials > 1
          ? IsLevelsetPeriodic(boundary_condition_reader, Direction::X)
          : false;
  bool material_periodic =
      IsMaterialPeriodic(boundary_condition_reader, Direction::X);

  if (levelset_periodic || material_periodic) {
    if (levelset_periodic != material_periodic && number_of_materials > 1) {
      throw std::invalid_argument(
          "Incorrect use of East-West periodic condition, both boundaries from "
          "the levelset and material must be periodic!");
    } else {
      active_periodic_locations |= PeriodicBoundariesLocations::EastWest;
    }
  }

  // Check the y-direction (levelset only if multi material)
  if constexpr (CC::DIM() != Dimension::One) {
    levelset_periodic =
        number_of_materials > 1
            ? IsLevelsetPeriodic(boundary_condition_reader, Direction::Y)
            : false;
    material_periodic =
        IsMaterialPeriodic(boundary_condition_reader, Direction::Y);

    if (levelset_periodic || material_periodic) {
      if (levelset_periodic != material_periodic && number_of_materials > 1) {
        throw std::invalid_argument(
            "Incorrect use of North-South periodic condition, both boundaries "
            "from the levelset and material must be periodic!");
      } else {
        active_periodic_locations |= PeriodicBoundariesLocations::NorthSouth;
      }
    }
  }

  // Check the z-direction (levelset only if multi material)
  if constexpr (CC::DIM() == Dimension::Three) {
    levelset_periodic =
        number_of_materials > 1
            ? IsLevelsetPeriodic(boundary_condition_reader, Direction::Z)
            : false;
    material_periodic =
        IsMaterialPeriodic(boundary_condition_reader, Direction::Z);

    if (levelset_periodic || material_periodic) {
      if (levelset_periodic != material_periodic && number_of_materials > 1) {
        throw std::invalid_argument(
            "Incorrect use of Top-Bottom periodic condition, both boundaries "
            "from the levelset and material must be periodic!");
      } else {
        active_periodic_locations |= PeriodicBoundariesLocations::TopBottom;
      }
    }
  }

  return active_periodic_locations;
}

/**
 * @brief Instantiates the complete topology manager class with the given input
 * classes.
 * @param input_reader Reader that provides access to the full data of the input
 * file.
 * @param material_manager material_manager Instance providing instantiated
 * material data.
 * @return The fully instantiated TopologyManager class.
 */
TopologyManager
InstantiateTopologyManager(InputReader const &input_reader,
                           MaterialManager const &material_manager) {

  return TopologyManager(
      GetNumberOfNodesOnLevelZero(input_reader.GetMultiResolutionReader()),
      input_reader.GetMultiResolutionReader().ReadMaximumLevel(),
      GetActivePeriodicDirections(input_reader.GetBoundaryConditionReader(),
                                  material_manager));
}
} // namespace Instantiation
