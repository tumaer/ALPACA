//===------------------- interface_tag_functions.cpp ----------------------===//
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
#include "interface_tag_functions.h"

#include "enums/interface_tag_definition.h"
#include "levelset/geometry/geometry_calculator.h"
#include "user_specifications/two_phase_constants.h"
#include "utilities/mathematical_functions.h"

namespace InterfaceTagFunctions {

/**
 * @brief Initializes all tags by default as cut-cells. This is necessary for
 * the function InterfaceTagManager::SetInternalCutCellTagsFromLevelset.
 * @param interface_tags Indirect return parameter for the interface tags.
 */
void InitializeInternalInterfaceTags(
    std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]) {

  for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        interface_tags[i][j][k] = ITTI(IT::OldCutCell);
      } // k
    }   // j
  }     // i
}

/**
 * @brief Sets the internal cut-cell tags according to the internal cells of the
 * given levelset buffer.
 * @param levelset Reference to the levelset values.
 * @param interface_tags Indirect return parameter for the derived interface
 * tags.
 */
void SetInternalCutCellTagsFromLevelset(
    double const (&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()],
    std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]) {

  // find cut-cells in internal cells - due to the CFL condition, it is enough
  // to test only those cells which were cut cells or neighbors before of cut
  // cells
  for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        // cell was close to the interface in previous step
        if (std::abs(interface_tags[i][j][k]) < ITTI(IT::ExtensionBand) ||
            std::abs(interface_tags[i][j][k]) == ITTI(IT::ScaleSeparatedCell)) {
          if (IsCutCell<GeometryCalculationSettings::CutCellCriteria>(
                  levelset, i, j, k)) {
            // this cut-cell was a cut-cell in the previous step
            if (std::abs(interface_tags[i][j][k]) <= ITTI(IT::NewCutCell))
              interface_tags[i][j][k] = ITTI(IT::OldCutCell);
            // this cut-cell was not a cut-cell before, but is now
            else
              interface_tags[i][j][k] =
                  Signum(levelset[i][j][k]) * ITTI(IT::NewCutCell);
          } else {
            // this is not a cut-cell -> reset to bulk phase. Later, this cell
            // interface tag is correctly set in function
            // InterfaceTagManager::SetTotalInterfaceTagsFromCutCells.
            interface_tags[i][j][k] =
                Signum(levelset[i][j][k]) * ITTI(IT::BulkPhase);
          }
        } else {
          // this cannot be a cut-cell, as it was too far away from the
          // interface -> reset it to bulk phase. Later, this cell interface tag
          // is correctly set in function
          // InterfaceTagManager::SetTotalInterfaceTagsFromCutCells.
          interface_tags[i][j][k] =
              Signum(levelset[i][j][k]) * ITTI(IT::BulkPhase);
        }
      } // k
    }   // j
  }     // i
}

/**
 * @brief Updates the interface tags with respect to all cut cells within the
 * block.
 * @param interface_tags Reference to the interface tag buffer that has to be
 * updated and already holds the cut cell tags.
 */
void SetTotalInterfaceTagsFromCutCells(
    std::int8_t (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]) {

  // these type-casts are necessary to allow for checking for cut-cells in the
  // outermost halo cell
  constexpr int total_cells_x = static_cast<int>(CC::TCX());
  constexpr int total_cells_y = static_cast<int>(CC::TCY());
  constexpr int total_cells_z = static_cast<int>(CC::TCZ());
  constexpr int cut_cell_neighbour_band_width_x = 1;
  constexpr int cut_cell_neighbour_band_width_y =
      CC::DIM() != Dimension::One ? 1 : 0;
  constexpr int cut_cell_neighbour_band_width_z =
      CC::DIM() == Dimension::Three ? 1 : 0;
  constexpr int extension_band_width_x = static_cast<int>(CC::EBW());
  constexpr int extension_band_width_y =
      CC::DIM() != Dimension::One ? static_cast<int>(CC::EBW()) : 0;
  constexpr int extension_band_width_z =
      CC::DIM() == Dimension::Three ? static_cast<int>(CC::EBW()) : 0;
  constexpr int reinitialization_band_width_x = static_cast<int>(CC::RBW());
  constexpr int reinitialization_band_width_y =
      CC::DIM() != Dimension::One ? static_cast<int>(CC::RBW()) : 0;
  constexpr int reinitialization_band_width_z =
      CC::DIM() == Dimension::Three ? static_cast<int>(CC::RBW()) : 0;

  // set interface tags from cut-cells in full block
  for (int i = 0; i < total_cells_x; ++i) {
    for (int j = 0; j < total_cells_y; ++j) {
      for (int k = 0; k < total_cells_z; ++k) {

        // find all cells that are now cut cells
        if (std::abs(interface_tags[i][j][k]) <= ITTI(IT::NewCutCell)) {

          // tag cells which are neighbors to cut cells
          for (int r = std::max(i - cut_cell_neighbour_band_width_x, 0);
               r <=
               std::min(i + cut_cell_neighbour_band_width_x, total_cells_x - 1);
               ++r) {
            for (int s = std::max(j - cut_cell_neighbour_band_width_y, 0);
                 s <= std::min(j + cut_cell_neighbour_band_width_y,
                               total_cells_y - 1);
                 ++s) {
              for (int t = std::max(k - cut_cell_neighbour_band_width_z, 0);
                   t <= std::min(k + cut_cell_neighbour_band_width_z,
                                 total_cells_z - 1);
                   ++t) {
                if (std::abs(interface_tags[r][s][t]) > ITTI(IT::NewCutCell))
                  interface_tags[r][s][t] = Signum(interface_tags[r][s][t]) *
                                            ITTI(IT::CutCellNeighbor);
              } // t
            }   // s
          }     // r

          // tag the extension band
          for (int r = std::max(i - extension_band_width_x, 0);
               r <= std::min(i + extension_band_width_x, total_cells_x - 1);
               ++r) {
            for (int s = std::max(j - extension_band_width_y, 0);
                 s <= std::min(j + extension_band_width_y, total_cells_y - 1);
                 ++s) {
              for (int t = std::max(k - extension_band_width_z, 0);
                   t <= std::min(k + extension_band_width_z, total_cells_z - 1);
                   ++t) {
                if (std::abs(interface_tags[r][s][t]) >
                    ITTI(IT::CutCellNeighbor))
                  interface_tags[r][s][t] =
                      Signum(interface_tags[r][s][t]) * ITTI(IT::ExtensionBand);
              } // t
            }   // s
          }     // r

          // tag the reinitialization band
          for (int r = std::max(i - reinitialization_band_width_x, 0);
               r <=
               std::min(i + reinitialization_band_width_x, total_cells_x - 1);
               ++r) {
            for (int s = std::max(j - reinitialization_band_width_y, 0);
                 s <=
                 std::min(j + reinitialization_band_width_y, total_cells_y - 1);
                 ++s) {
              for (int t = std::max(k - reinitialization_band_width_z, 0);
                   t <= std::min(k + reinitialization_band_width_z,
                                 total_cells_z - 1);
                   ++t) {
                if (std::abs(interface_tags[r][s][t]) > ITTI(IT::ExtensionBand))
                  interface_tags[r][s][t] = Signum(interface_tags[r][s][t]) *
                                            ITTI(IT::ReinitializationBand);
              } // t
            }   // s
          }     // r
        }
      } // k
    }   // j
  }     // i
}

/**
 * @brief Gives whether or not the given interface tags are uniform (hence, a
 * bulk phase).
 * @param interface_tags Reference to the interface tags to be checked.
 * @return True if all interface tags are uniform, false otherwise.
 */
bool TotalInterfaceTagsAreUniform(
    std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()]) {
  for (unsigned int i = 0; i < CC::TCX(); ++i) {
    for (unsigned int j = 0; j < CC::TCY(); ++j) {
      for (unsigned int k = 0; k < CC::TCZ(); ++k) {
        if (std::abs(interface_tags[i][j][k]) != ITTI(IT::BulkPhase)) {
          return false;
        }
      } // k
    }   // j
  }     // i
  return true;
}

} // namespace InterfaceTagFunctions
