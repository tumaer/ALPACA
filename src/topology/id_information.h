//===------------------------- id_information.h ---------------------------===//
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
#ifndef ID_INFORMATION_H
#define ID_INFORMATION_H

#include "boundary_condition/boundary_specifications.h"
#include "topology/node_id_type.h"
#include "user_specifications/compile_time_constants.h"
#include <array>
#include <vector>

nid_t EastNeighborOfNodeWithId(nid_t const id);
nid_t WestNeighborOfNodeWithId(nid_t const id);
nid_t NorthNeighborOfNodeWithId(nid_t const id);
nid_t SouthNeighborOfNodeWithId(nid_t const id);
nid_t TopNeighborOfNodeWithId(nid_t const id);
nid_t BottomNeighborOfNodeWithId(nid_t const id);

nid_t GetNeighborId(nid_t const id, BoundaryLocation const location);

nid_t IdSeed();
nid_t CutHeadBit(nid_t id, unsigned int const level);
nid_t AddHeadBit(nid_t const id, unsigned int const level);

bool IsNaturalExternalBoundary(
    BoundaryLocation const location, nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero);
bool IsExternalBoundary(
    BoundaryLocation const location, nid_t const id,
    std::array<unsigned int, 3> const number_of_nodes_on_level_zero);

nid_t RawMortonIndexOfId(nid_t const id);

/**
 * @brief Gives the level that the given id is associated to. $NOT SAFE, WRONG
 * INPUT RESULTS IN QUITE INCORRECT DATA$.
 * @param id Id of node to be evaluated.
 * @return Level of the id.
 */
constexpr unsigned int LevelOfNode(nid_t const id) {

  nid_t working_id = (id >> 21); // Cut Shadows
  unsigned int level = 0;

  while (working_id > 10) {
    working_id >>= 3;
    level++;
  }

  return level;
}

/**
 * @brief Gives the cell size for a node a given level.
 * @param node_size_on_level_zero The size (= size of internal cells) of a node
 * on level zero.
 * @param level The level for which the cell size is evaluated.
 * @return The cell size of a node on the given level.
 */
constexpr double CellSizeOfLevel(double const node_size_on_level_zero,
                                 unsigned int const level) {
  return node_size_on_level_zero / double(CC::ICX()) / double(1 << level);
}

/**
 * @brief Gives the cell size for the node of the given id.
 * @param id Id of the node whose cell size is to be evaluated.
 * @param node_size_on_level_zero The size (= size of internal cells) of a node
 * on level zero.
 * @return The cell size of the node.
 */
constexpr double CellSizeOfId(nid_t const id,
                              double const node_size_on_level_zero) {
  return CellSizeOfLevel(node_size_on_level_zero, LevelOfNode(id));
}

/**
 * @brief Gives the length of the internal domain within the node of the given
 * id.
 * @param id Id of the node whose size is to be evaluated.
 * @param node_size_on_level_zero The size (= size of internal cells) of a node
 * on level zero.
 * @return Length of the internal domain in the node.
 */
constexpr double DomainSizeOfId(nid_t const id,
                                double const node_size_on_level_zero) {
  nid_t divisor = 1 << (LevelOfNode(id));
  return node_size_on_level_zero / double(divisor);
}

/**
 * @brief Gives the coordinates of the coordinate system origin of a nodes
 * internal cells, i.e. cell [CC::FICX()][CC::FICY()][CC::FICZ()].
 * @param id The id of a node for which the coordinates are calculated.
 * @param block_size The size (= size of internal cells) of a block.
 * @return The coordinates of the coordinate origin in the block, i.e.
 * coordinates of corner of first internal cells.
 */
constexpr std::array<double, 3> DomainCoordinatesOfId(nid_t const id,
                                                      double const block_size) {
  nid_t level_operator = id;
  unsigned int last_bit = (level_operator & 0x1);
  double size = block_size;

  double x = 0;
  double y = 0;
  double z = 0;

  while (level_operator != 10) {
    x += last_bit * size;
    level_operator >>= 1; // Shift by one to the right
    last_bit = (level_operator & 0x1);
    y += last_bit * size;
    level_operator >>= 1; // Shift by one to the right
    last_bit = (level_operator & 0x1);
    z += last_bit * size;
    level_operator >>= 1; // Shift by one to the right
    last_bit = (level_operator & 0x1);
    size *= 2;
  }

  return {x, y, z};
}

/**
 * @brief Gives the unique identifer of the parent node. $NO INPUT CHECKS,
 * CALLER MUST ENSURE CORRECT INPUT$.
 * @param id Id of the Child.
 * @return Id of the Parent.
 */
constexpr nid_t ParentIdOfNode(nid_t const id) { return (id >> 3); }

/**
 * @brief Indicates the position of the Node among its siblings. The position is
 * encoded as integer value with 0 = bottom-south-west to 7 = top-north-east.
 * @param id The id of the node, whose position is to be determined.
 * @return The position among its siblings.
 */
constexpr unsigned int PositionOfNodeAmongSiblings(nid_t const id) {
  return (id & 0x7);
}

/**
 * @brief Functions to indicate whether or not a node is (one of the) east
 * (other locations respectively) most among its siblings. I.e. has (at least
 * one) a sibling to its west (other locations respectively).
 * @param id The id ot the node.
 * @return True if it is most east (other locations respectively), False
 * otherwise.
 */
constexpr bool EastInSiblingPack(nid_t const id) { return (id & 0x1) == 1; }
constexpr bool WestInSiblingPack(nid_t const id) { return (id & 0x1) == 0; }
constexpr bool NorthInSiblingPack(nid_t const id) { return (id & 0x2) == 2; }
constexpr bool SouthInSiblingPack(nid_t const id) { return (id & 0x2) == 0; }
constexpr bool TopInSiblingPack(nid_t const id) { return (id & 0x4) == 4; }
constexpr bool BottomInSiblingPack(nid_t const id) { return (id & 0x4) == 0; }
/**@}*/

/**
 * @brief Gives the Ids of all eight children of a Node.
 * @param id The Id of the parent node
 * @return Ids of the children in increasing order, i.e. bottom-south-west,
 * bottom-south-east, bottom-north-west, ... ,top-north-east.
 */
inline std::vector<nid_t> IdsOfChildren(nid_t const parent_id) {
  return {{(parent_id << 3), (parent_id << 3) + 1
#if DIMENSION > 1
           ,
           (parent_id << 3) + 2, (parent_id << 3) + 3
#endif
#if DIMENSION == 3
           ,
           (parent_id << 3) + 4, (parent_id << 3) + 5, (parent_id << 3) + 6,
           (parent_id << 3) + 7
#endif
  }};
}

#endif // ID_INFORMATION_H
