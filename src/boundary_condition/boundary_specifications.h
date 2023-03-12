//===-------------------- boundary_specifications.h -----------------------===//
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
#ifndef BOUNDARY_SPECIFICATIONS_H
#define BOUNDARY_SPECIFICATIONS_H

#include <stdexcept>
#include <string>
#include <type_traits>

#include "utilities/string_operations.h"

/**
 * @brief Unique identifier to separate the (implemented) types of material
 * boundary condition.
 */
enum class MaterialBoundaryType {
  Internal,
  ZeroGradient,
  Symmetry,
  FixedValue,
  Wall,
  Periodic
};

/**
 * @brief Unique identifier to separate the (implemented) types of level-set
 * boundary condition.
 */
enum class LevelSetBoundaryType { Internal, ZeroGradient, Symmetry, Periodic };

/**
 * @brief Unique Identifier for the Boundary Location (Compass orientations plus
 * Top/Bottom plus Diagonals)
 * @note Do NOT change the underlying type and index structure. This would break
 * the algorithm in multiple places!
 */
enum class BoundaryLocation : unsigned int {
  // normal Types used in 1D, for planes, and for DomainBoundaries
  East = 0,
  West = 1,
  North = 2,
  South = 3,
  Top = 4,
  Bottom = 5,
  // Diagonals in 2 Dimensions for Sticks in 3D and cubes in 2D
  BottomNorth,
  BottomSouth,
  TopNorth,
  TopSouth, // x-Axis Sticks
  BottomEast,
  BottomWest,
  TopEast,
  TopWest, // y-Axis Sticks
  NorthEast,
  NorthWest,
  SouthEast,
  SouthWest, // z-Axis Sticks
             // Diagonals in 3 Dimensions for cubes in 3D
  EastNorthTop,
  EastNorthBottom,
  EastSouthTop,
  EastSouthBottom, // East-Side
  WestNorthTop,
  WestNorthBottom,
  WestSouthTop,
  WestSouthBottom, // West-Side
  // DO NOT EXTEND WITHOUT ADJUSTING THE OppositeDirection FUNCTION
};

/**
 * @brief Converts a BoundaryLocation identifier to a (C++11 standard compliant,
 * i. e. positive) array index. "LTI = Location To Index".
 * @param l The location identifier.
 * @return Index to be used in Arrays.
 */
static inline constexpr std::underlying_type<BoundaryLocation>::type
LTI(BoundaryLocation const l) {
  return static_cast<typename std::underlying_type<BoundaryLocation>::type>(l);
}

/**
 * @brief Gives the inverse direction. Convenience function.
 * @param location The direction that is to be inverted.
 * @return The inverse direction.
 */
constexpr BoundaryLocation OppositeDirection(BoundaryLocation const location) {
  switch (location) {
  case BoundaryLocation::East: // Planes
    return BoundaryLocation::West;
  case BoundaryLocation::West:
    return BoundaryLocation::East;
  case BoundaryLocation::North:
    return BoundaryLocation::South;
  case BoundaryLocation::South:
    return BoundaryLocation::North;
  case BoundaryLocation::Top:
    return BoundaryLocation::Bottom;
  case BoundaryLocation::Bottom:
    return BoundaryLocation::Top;

  case BoundaryLocation::BottomNorth: // Sticks
    return BoundaryLocation::TopSouth;
  case BoundaryLocation::BottomSouth:
    return BoundaryLocation::TopNorth;
  case BoundaryLocation::TopNorth:
    return BoundaryLocation::BottomSouth;
  case BoundaryLocation::TopSouth:
    return BoundaryLocation::BottomNorth;

  case BoundaryLocation::BottomEast:
    return BoundaryLocation::TopWest;
  case BoundaryLocation::BottomWest:
    return BoundaryLocation::TopEast;
  case BoundaryLocation::TopEast:
    return BoundaryLocation::BottomWest;
  case BoundaryLocation::TopWest:
    return BoundaryLocation::BottomEast;

  case BoundaryLocation::NorthEast:
    return BoundaryLocation::SouthWest;
  case BoundaryLocation::NorthWest:
    return BoundaryLocation::SouthEast;
  case BoundaryLocation::SouthEast:
    return BoundaryLocation::NorthWest;
  case BoundaryLocation::SouthWest:
    return BoundaryLocation::NorthEast;

  case BoundaryLocation::EastNorthTop: // Cubes
    return BoundaryLocation::WestSouthBottom;
  case BoundaryLocation::EastNorthBottom:
    return BoundaryLocation::WestSouthTop;
  case BoundaryLocation::EastSouthTop:
    return BoundaryLocation::WestNorthBottom;
  case BoundaryLocation::EastSouthBottom:
    return BoundaryLocation::WestNorthTop;
  case BoundaryLocation::WestNorthTop:
    return BoundaryLocation::EastSouthBottom;
  case BoundaryLocation::WestNorthBottom:
    return BoundaryLocation::EastSouthTop;
  case BoundaryLocation::WestSouthTop:
    return BoundaryLocation::EastNorthBottom;
  default: // Only remaning case: BoundaryLocation::WestSouthBottom:
    return BoundaryLocation::EastNorthTop;
  }
}

/**
 * @brief Gives the string for each boundary location.
 * @param location Boundary location for which the string should be returned.
 * @param first_capitalized Flag whether the first letter should be capitalized
 * or not.
 * @return string of boundary location.
 *
 * @note Do not change leading Upper case letter. Required for rading of
 * boundary conditions.
 */
inline std::string BoundaryLocationToString(BoundaryLocation const location,
                                            bool const first_capitalized) {
  switch (location) {
  case BoundaryLocation::East: {
    return first_capitalized ? "East" : "east";
  }
  case BoundaryLocation::West: {
    return first_capitalized ? "West" : "west";
  }
  case BoundaryLocation::North: {
    return first_capitalized ? "North" : "north";
  }
  case BoundaryLocation::South: {
    return first_capitalized ? "South" : "south";
  }
  case BoundaryLocation::Top: {
    return first_capitalized ? "Top" : "top";
  }
  case BoundaryLocation::Bottom: {
    return first_capitalized ? "Bottom" : "bottom";
  }
  default: {
    throw std::logic_error("Boundary location not known!");
  }
  }
}

/**
 * @brief Gives the proper Material boundary type for a given string.
 * @param boundary_string String that should be converted.
 * @return Material boundary type identifier.
 */
inline MaterialBoundaryType
StringToMaterialBoundaryType(std::string const &boundary_string) {
  // transform string to upper case without spaces
  std::string const boundary_upper_case(
      StringOperations::ToUpperCaseWithoutSpaces(boundary_string));
  // switch statements cannot be used with strings
  if (boundary_upper_case == "ZEROGRADIENT") {
    return MaterialBoundaryType::ZeroGradient;
  } else if (boundary_upper_case == "SYMMETRY") {
    return MaterialBoundaryType::Symmetry;
  } else if (boundary_upper_case == "WALL") {
    return MaterialBoundaryType::Wall;
  } else if (boundary_upper_case == "FIXEDVALUE") {
    return MaterialBoundaryType::FixedValue;
  } else if (boundary_upper_case == "PERIODIC") {
    return MaterialBoundaryType::Periodic;
  } else {
    throw std::logic_error("Material boundary type '" + boundary_upper_case +
                           "' not known!");
  }
}

/**
 * @brief Gives the proper Levelset boundary type for a given string
 * @param boundary_string String that should be converted
 * @return Levelset boundary type identifier
 */
inline LevelSetBoundaryType
StringToLevelSetBoundaryType(std::string const &boundary_string) {
  // transform string to upper case without spaces
  std::string const boundary_upper_case(
      StringOperations::ToUpperCaseWithoutSpaces(boundary_string));
  // switch statements cannot be used with strings
  if (boundary_upper_case == "ZEROGRADIENT") {
    return LevelSetBoundaryType::ZeroGradient;
  } else if (boundary_upper_case == "SYMMETRY") {
    return LevelSetBoundaryType::Symmetry;
  } else if (boundary_upper_case == "PERIODIC") {
    return LevelSetBoundaryType::Periodic;
  } else {
    throw std::logic_error("Levelset boundary type '" + boundary_upper_case +
                           "' not known!");
  }
}

#endif // BOUNDARY_SPECIFICATIONS_H
