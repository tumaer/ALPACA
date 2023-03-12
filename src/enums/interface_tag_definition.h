//===-------------------- interface_tag_definition.h ----------------------===//
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
#ifndef INTERFACE_TAG_DEFINITION_H
#define INTERFACE_TAG_DEFINITION_H

#include <cstdint>
#include <type_traits>

/**
 * This info is to stay during development in case we need to further update the
 * tagging-system, as currently functionality can only partially be tested. Old
 * Cut Cell           = 0 -> already cut-cell in last RK-step, still cut-cell.
 * Is integrated in both RK steps. New Cut Cell           = 1 -> no cut-cell in
 * last RK-step, now a cut-cell. TODO-19 Check whether they can be neglected for
 * integration in second RK step. Would be filled by mixing. = 2 -> was a
 * cut-cell in last RK-step, now empty/full cell close to interface TODO-19 JK
 * not yet implemented. Could be helpful for mixing to avoid violation of
 * conservation once cells should be empty but still contain conservatives.
 *
 * Cut-cell neighbour     = 3  -> indicates whether cell has a cut-cell
 * neighbour, tags 0, 1 and 2 are dominating.
 *
 *                        = 4  -> was a cut-cell in last RK-step, now empty/full
 * cell away from interface.
 *
 * Extension-band         = 7  -> indicates whether material is extended in this
 * region, tags 0, 1, 2 and 3 are dominating, bit-shift gives width of extension
 * (3 currently), keep this in mind when changing extension width!
 *
 * Reinitialization-band  = 8  -> indicates whether levelset is reinitialized in
 * this region, usually one cell-layer more than extension-width, depends on the
 * order of the normal-calculation finite-difference scheme. Bit-shift gives
 * width of reinitialization (4 currently), keep this in mind when changing
 * reinitialization width!
 *
 * Bulk-face              = 10 -> indicates a bulk-face cell, tags 0, 1, 2, 3, 7
 * and 8 are dominating.
 *
 * Scale-separation       = 50 -> whether a cell is scale-separated, only
 * used/set in scale-separation, overwritten afterwards in SetInterfaceTags.
 */

/**
 * @brief Identifiers for the interface tags to be used.
 *
 * @note It is absolutely neccesary to keep cut cell < 2, extension band <= 7,
 * Reinitialization band <= 8, and bulk phase and scale separated cells > 10 for
 * operations acting on these enums. Values for extension and reinitialization
 * are used as they give the currently used values for extension and
 * reinitialization (3 and 4) by a bit-shift to the right. If this needs to be
 * changed, it can be done here.
 */
enum class InterfaceTag : std::int8_t {
  OldCutCell = 0, // this cell was a cut-cell in the last iteration and still is
  NewCutCell = 1, // this cell is now a cut-cell, but was not before
  CutCellNeighbor = 3, // this cell has a cut-cell neighbor, i.e. its levelset
                       // is advected but not reinitialized
  ExtensionBand =
      7, // the ghost material of these cells needs to be filled by extension
  ReinitializationBand =
      8,          // the level set of these cells needs to be reinitialized
  BulkPhase = 10, // these are cells far away from the interface
  ScaleSeparatedCell = 50
}; // these cells were changed during the scale-separation algorithm

/**
 * @brief Converts an interface-tag identifier to a (C++11 standard compliant,
 * i. e. positive) array index. "ITTI = Interface tag To Index"
 * @param it The interface-tag identifier.
 * @return Integer to be used in arrays.
 */
constexpr std::underlying_type<InterfaceTag>::type ITTI(InterfaceTag const it) {
  return static_cast<typename std::underlying_type<InterfaceTag>::type>(it);
}

using IT = InterfaceTag;

#endif // INTERFACE_TAG_DEFINITION_H
