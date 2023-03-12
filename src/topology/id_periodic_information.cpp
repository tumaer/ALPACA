//===----------------- id_periodic_information.cpp ------------------------===//
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
#include "topology/id_periodic_information.h"
#include "topology/id_information.h"
#include <bitset>

namespace {

/**
 * @brief Checks whether a given neighbor location of a node represents an
 * external periodic location.
 * @param periodic_location Periodic location to be checked.
 * @param neighbor_location Location where the neighbor node lies.
 * @param id Node if that is investigated.
 * @param number_of_nodes_on_level_zero Number of nodes on level zero in the
 * three Cartesian directions.
 * @param active_periodic_locations Active periodic locations for the simulation
 * (1: East-West, 2:North-South, 4:Top-Bottom).
 * @return True if neighbor location is an external periodic boundary, otherwise
 * False.
 */
bool NeighborIsExternalPeriodic(
    PeriodicBoundariesLocations const periodic_location,
    BoundaryLocation const neighbor_location, nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {
  return (active_periodic_locations & periodic_location) &&
         IsNaturalExternalBoundary(neighbor_location, id,
                                   number_of_nodes_on_level_zero);
}

/**
 * @brief Calculates the id of the eastern periodic neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @param active_periodic_locations bitwise representation of the active
 * periodic locations.
 * @return Id of eastern periodic neighbor.
 */
nid_t EastPeriodicNeighborOfNodeWithId(
    nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {

  if (NeighborIsExternalPeriodic(
          PeriodicBoundariesLocations::EastWest, BoundaryLocation::East, id,
          number_of_nodes_on_level_zero, active_periodic_locations)) {
    unsigned int const level = LevelOfNode(id);
    std::bitset<64> const east_mask(0xDB6DB6DB6DB6DB6);
    std::bitset<64> working_id(CutHeadBit(id, level));
    working_id &= east_mask;
    return (AddHeadBit(working_id.to_ullong(), level));
  } else {
    return EastNeighborOfNodeWithId(id);
  }
}

/**
 * @brief Calculates the id of the western periodic neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @param active_periodic_locations bitwise representation of the active
 * periodic locations.
 * @return Id of western periodic neighbor.
 */
nid_t WestPeriodicNeighborOfNodeWithId(
    nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {

  if (NeighborIsExternalPeriodic(
          PeriodicBoundariesLocations::EastWest, BoundaryLocation::West, id,
          number_of_nodes_on_level_zero, active_periodic_locations)) {
    unsigned int const level = LevelOfNode(id);
    std::bitset<64> working_id(CutHeadBit(id, level));
    std::bitset<64> const x_coordinate(
        number_of_nodes_on_level_zero[0] * (1u << level) - 1u);
    for (std::uint8_t i = 0; i < 60; i += 3) {
      working_id[i] = x_coordinate[i / 3];
    }
    return (AddHeadBit(working_id.to_ullong(), level));
  } else {
    return WestNeighborOfNodeWithId(id);
  }
}

/**
 * @brief Calculates the id of the northern periodic neighbor of the id
 * provided.
 * @param id The id of the node whose neighbor is to be found.
 * @param active_periodic_locations bitwise representation of the active
 * periodic locations.
 * @return Id of northern periodic neighbor.
 */
nid_t NorthPeriodicNeighborOfNodeWithId(
    nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {

  if (NeighborIsExternalPeriodic(
          PeriodicBoundariesLocations::NorthSouth, BoundaryLocation::North, id,
          number_of_nodes_on_level_zero, active_periodic_locations)) {
    unsigned int const level = LevelOfNode(id);
    std::bitset<64> const north_mask(0x5B6DB6DB6DB6DB6D);
    std::bitset<64> working_id(CutHeadBit(id, level));
    working_id &= north_mask;
    return (AddHeadBit(working_id.to_ullong(), level));
  } else {
    return NorthNeighborOfNodeWithId(id);
  }
}

/**
 * @brief Calculates the id of the southern  periodic neighbor of the id
 * provided.
 * @param id The id of the node whose neighbor is to be found.
 * @param active_periodic_locations bitwise representation of the active
 * periodic locations.
 * @return Id of southern  periodic neighbor.
 */
nid_t SouthPeriodicNeighborOfNodeWithId(
    nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {

  if (NeighborIsExternalPeriodic(
          PeriodicBoundariesLocations::NorthSouth, BoundaryLocation::South, id,
          number_of_nodes_on_level_zero, active_periodic_locations)) {
    unsigned int const level = LevelOfNode(id);
    std::bitset<64> const y_coordinate(
        number_of_nodes_on_level_zero[1] * (1u << level) - 1u);
    std::bitset<64> working_id(CutHeadBit(id, level));
    for (std::uint8_t i = 1; i < 60; i += 3) {
      working_id[i] = y_coordinate[i / 3];
    }
    return (AddHeadBit(working_id.to_ullong(), level));
  } else {
    return SouthNeighborOfNodeWithId(id);
  }
}

/**
 * @brief Calculates the id of the top periodic neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @param active_periodic_locations bitwise representation of the active
 * periodic locations.
 * @return Id of top periodic neighbor.
 */
nid_t TopPeriodicNeighborOfNodeWithId(
    nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {

  if (NeighborIsExternalPeriodic(
          PeriodicBoundariesLocations::TopBottom, BoundaryLocation::Top, id,
          number_of_nodes_on_level_zero, active_periodic_locations)) {
    unsigned int const level = LevelOfNode(id);
    std::bitset<64> const top_mask(0x36DB6DB6DB6DB6DB);
    std::bitset<64> working_id(CutHeadBit(id, level));
    working_id &= top_mask;
    return (AddHeadBit(working_id.to_ullong(), level));
  } else {
    return (TopNeighborOfNodeWithId(id));
  }
}

/**
 * @brief Calculates the id of the bottom periodic neighbor of the id provided.
 * @param id The id of the node whose neighbor is to be found.
 * @param active_periodic_locations bitwise representation of the active
 * periodic locations.
 * @return Id of bottom periodic neighbor.
 */
nid_t BottomPeriodicNeighborOfNodeWithId(
    nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {

  if (NeighborIsExternalPeriodic(
          PeriodicBoundariesLocations::TopBottom, BoundaryLocation::Bottom, id,
          number_of_nodes_on_level_zero, active_periodic_locations)) {
    unsigned int const level = LevelOfNode(id);
    std::bitset<64> working_id(CutHeadBit(id, level));
    std::bitset<64> const z_coordinate(
        number_of_nodes_on_level_zero[2] * (1u << level) - 1u);
    for (std::uint8_t i = 2; i < 60; i += 3) {
      working_id[i] = z_coordinate[i / 3];
    }
    return (AddHeadBit(working_id.to_ullong(), level));
  } else {
    return (BottomNeighborOfNodeWithId(id));
  }
}

/**
 * @brief Determines whether a node face is also an external or domain edge.
 * @param location The  direction of the edge under consideration.
 * @param id The id of the node under investigation.
 * @param setup The simulation settings as provided by the user.
 * @return True if the edge is a domain edge, false otherwise, i.e. internal
 * edge.
 */
bool PeriodicIsNaturalExternalBoundary(
    BoundaryLocation const location, nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {

  std::bitset<64> const input(CutHeadBit(id, LevelOfNode(id)));

  switch (location) {
  case BoundaryLocation::West: {
    std::bitset<64> const west_mask(0x249249249249249);
    std::bitset<64> const result(input & west_mask);
    if (!result.to_ullong() &&
        !(active_periodic_locations & PeriodicBoundariesLocations::EastWest)) {
      return true;
    }
  } break;
  case BoundaryLocation::South: {
    std::bitset<64> const south_mask(0x492492492492492);
    std::bitset<64> const result(input & south_mask);
    if (!result.to_ullong() && !(active_periodic_locations &
                                 PeriodicBoundariesLocations::NorthSouth)) {
      return true;
    }
  } break;
  case BoundaryLocation::Bottom: {
    std::bitset<64> const bottom_mask(0x924924924924924);
    std::bitset<64> const result(input & bottom_mask);
    if (!result.to_ullong() &&
        !(active_periodic_locations & PeriodicBoundariesLocations::TopBottom)) {
      return true;
    }
  } break;
  case BoundaryLocation::East: {
    std::uint32_t const tows_exponent(1u << LevelOfNode(id));
    std::bitset<20> const east_mask(
        (number_of_nodes_on_level_zero[0] * tows_exponent) - 1u);
    std::bitset<1> indicator(1);
    for (unsigned int i = 0; i < 20; ++i) {
      indicator[0] = (~(input[i * 3] ^ east_mask[i]) & indicator[0]);
    }
    if (indicator[0] &&
        !(active_periodic_locations & PeriodicBoundariesLocations::EastWest)) {
      return true;
    }
  } break;
  case BoundaryLocation::North: {
    std::uint32_t const tows_exponent(1u << LevelOfNode(id));
    std::bitset<20> const north_mask(
        (number_of_nodes_on_level_zero[1] * tows_exponent) - 1u);
    std::bitset<1> indicator(1);
    for (unsigned int i = 0; i < 20; ++i) {
      indicator[0] = (~(input[(i * 3) + 1] ^ north_mask[i]) & indicator[0]);
    }
    if (indicator[0] && !(active_periodic_locations &
                          PeriodicBoundariesLocations::NorthSouth)) {
      return true;
    }
  } break;
#ifndef PERFORMANCE
  case BoundaryLocation::Top: {
    std::uint32_t const tows_exponent(1u << LevelOfNode(id));
    std::bitset<20> top_mask(
        (number_of_nodes_on_level_zero[2] * tows_exponent) - 1u);
    std::bitset<1> indicator(1);
    for (unsigned int i = 0; i < 20; ++i) {
      indicator[0] = (~(input[(i * 3) + 2] ^ top_mask[i]) & indicator[0]);
    }
    if (indicator[0] &&
        !(active_periodic_locations & PeriodicBoundariesLocations::TopBottom)) {
      return true;
    }
  } break;
  default: {
    throw std::invalid_argument("Boundary Type in IsExternal does not exist");
  } break;
#else
  default: /* BoundaryLocation::Top */ {
    std::uint32_t const tows_exponent(1u << LevelOfNode(id));
    std::bitset<20> top_mask(
        (number_of_nodes_on_level_zero[2] * tows_exponent) - 1u);
    std::bitset<1> indicator(1);
    for (unsigned int i = 0; i < 20; ++i) {
      indicator[0] = (~(input[(i * 3) + 2] ^ top_mask[i]) & indicator[0]);
    }
    if (indicator[0] &&
        !(active_periodic_locations & PeriodicBoundariesLocations::TopBottom)) {
      return true;
    }
  }
#endif
  }

  return false;
}
} // namespace

/**
 * @brief Gives the id of a periodic neighbor at the provided direction.
 * @param id The id of the node whose neighbor is to be found.
 * @param location Direction in which the neighbor is located.
 * @param setup The simulation settings as provided by the user.
 * @param active_periodic_locations bitwise representation of the active
 * periodic locations.
 * @return Id of the periodic neighbor.
 */
nid_t GetPeriodicNeighborId(
    nid_t const id, BoundaryLocation const location,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {

  switch (location) {
  // Natural
  case BoundaryLocation::East:
    return EastPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                            active_periodic_locations);
  case BoundaryLocation::West:
    return WestPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                            active_periodic_locations);
  case BoundaryLocation::North:
    return NorthPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                             active_periodic_locations);
  case BoundaryLocation::South:
    return SouthPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                             active_periodic_locations);
  case BoundaryLocation::Top:
    return TopPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                           active_periodic_locations);
  case BoundaryLocation::Bottom:
    return BottomPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                              active_periodic_locations);

    // Sticks
  case BoundaryLocation::BottomNorth:
    return BottomPeriodicNeighborOfNodeWithId(
        NorthPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                          active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::BottomSouth:
    return BottomPeriodicNeighborOfNodeWithId(
        SouthPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                          active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::TopNorth:
    return TopPeriodicNeighborOfNodeWithId(
        NorthPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                          active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::TopSouth:
    return TopPeriodicNeighborOfNodeWithId(
        SouthPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                          active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);

  case BoundaryLocation::BottomEast:
    return BottomPeriodicNeighborOfNodeWithId(
        EastPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                         active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::BottomWest:
    return BottomPeriodicNeighborOfNodeWithId(
        WestPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                         active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::TopEast:
    return TopPeriodicNeighborOfNodeWithId(
        EastPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                         active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::TopWest:
    return TopPeriodicNeighborOfNodeWithId(
        WestPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                         active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);

  case BoundaryLocation::NorthEast:
    return NorthPeriodicNeighborOfNodeWithId(
        EastPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                         active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::NorthWest:
    return NorthPeriodicNeighborOfNodeWithId(
        WestPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                         active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::SouthEast:
    return SouthPeriodicNeighborOfNodeWithId(
        EastPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                         active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::SouthWest:
    return SouthPeriodicNeighborOfNodeWithId(
        WestPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                         active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);

    // Cubes
  case BoundaryLocation::EastNorthTop:
    return EastPeriodicNeighborOfNodeWithId(
        NorthPeriodicNeighborOfNodeWithId(
            TopPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                            active_periodic_locations),
            number_of_nodes_on_level_zero, active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::EastNorthBottom:
    return EastPeriodicNeighborOfNodeWithId(
        NorthPeriodicNeighborOfNodeWithId(
            BottomPeriodicNeighborOfNodeWithId(
                id, number_of_nodes_on_level_zero, active_periodic_locations),
            number_of_nodes_on_level_zero, active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::EastSouthTop:
    return EastPeriodicNeighborOfNodeWithId(
        SouthPeriodicNeighborOfNodeWithId(
            TopPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                            active_periodic_locations),
            number_of_nodes_on_level_zero, active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::EastSouthBottom:
    return EastPeriodicNeighborOfNodeWithId(
        SouthPeriodicNeighborOfNodeWithId(
            BottomPeriodicNeighborOfNodeWithId(
                id, number_of_nodes_on_level_zero, active_periodic_locations),
            number_of_nodes_on_level_zero, active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);

  case BoundaryLocation::WestNorthTop:
    return WestPeriodicNeighborOfNodeWithId(
        NorthPeriodicNeighborOfNodeWithId(
            TopPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                            active_periodic_locations),
            number_of_nodes_on_level_zero, active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::WestNorthBottom:
    return WestPeriodicNeighborOfNodeWithId(
        NorthPeriodicNeighborOfNodeWithId(
            BottomPeriodicNeighborOfNodeWithId(
                id, number_of_nodes_on_level_zero, active_periodic_locations),
            number_of_nodes_on_level_zero, active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::WestSouthTop:
    return WestPeriodicNeighborOfNodeWithId(
        SouthPeriodicNeighborOfNodeWithId(
            TopPeriodicNeighborOfNodeWithId(id, number_of_nodes_on_level_zero,
                                            active_periodic_locations),
            number_of_nodes_on_level_zero, active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  case BoundaryLocation::WestSouthBottom:
    return WestPeriodicNeighborOfNodeWithId(
        SouthPeriodicNeighborOfNodeWithId(
            BottomPeriodicNeighborOfNodeWithId(
                id, number_of_nodes_on_level_zero, active_periodic_locations),
            number_of_nodes_on_level_zero, active_periodic_locations),
        number_of_nodes_on_level_zero, active_periodic_locations);
  default:
    throw std::invalid_argument(
        "Boundary Location not Found in Node::GetPeriodicNeighborId() - "
        "Impossible Error");
  }
}

/**
 * @brief Wrapper for PeriodicIsNaturalExternalBoundary if natural is not
 * guaranteed. Determines whether a block edge is also an external or domain
 * edge.
 * @param location The  direction of the edge under consideration.
 * @param id The id of the node under investigation.
 * @param setup The simulation settings as provided by the user.
 * @return True if the edge is a domain edge, false otherwise, i.e. internal
 * edge.
 */
bool PeriodicIsExternalBoundary(
    BoundaryLocation const location, nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero,
    unsigned int const active_periodic_locations) {
  // natural | NH Such comparison are okay by (enforced) definiton of
  // BoundaryLocation
  if (LTI(location) <= LTI(BoundaryLocation::Bottom)) {
    return PeriodicIsNaturalExternalBoundary(
        location, id, number_of_nodes_on_level_zero, active_periodic_locations);
  }

  switch (location) {
  // Sticks
  case BoundaryLocation::BottomNorth:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::Bottom, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::North, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::BottomSouth:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::Bottom, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::South, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::TopNorth:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::Top, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::North, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::TopSouth:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::Top, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::South, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::BottomEast:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::Bottom, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::East, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::BottomWest:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::Bottom, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::West, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::TopEast:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::Top, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::East, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::TopWest:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::Top, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::West, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));

  case BoundaryLocation::NorthEast:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::North, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::East, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::NorthWest:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::North, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::West, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::SouthEast:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::South, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::East, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::SouthWest:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::South, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::West, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));

  // Cubes
  case BoundaryLocation::EastNorthTop:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::East, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::North, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::Top, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::EastNorthBottom:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::East, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::North, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::Bottom, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::EastSouthTop:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::East, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::South, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::Top, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::EastSouthBottom:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::East, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::South, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::Bottom, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::WestNorthTop:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::West, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::North, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::Top, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::WestNorthBottom:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::West, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::North, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::Bottom, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  case BoundaryLocation::WestSouthTop:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::West, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::South, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::Top, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
#ifndef PERFORMANCE
  case BoundaryLocation::WestSouthBottom:
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::West, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::South, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::Bottom, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
  default:
    throw std::invalid_argument(
        "Boundary Location not Found in GetNeighborId() - Impossible Error");
#else
  default: /* BoundaryLocation::WestSouthBottom */
    return (PeriodicIsNaturalExternalBoundary(BoundaryLocation::West, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::South, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations) ||
            PeriodicIsNaturalExternalBoundary(BoundaryLocation::Bottom, id,
                                              number_of_nodes_on_level_zero,
                                              active_periodic_locations));
#endif
  }
}
