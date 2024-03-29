//===------------------- aperture_cut_cell_mixer.cpp ----------------------===//
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
#include "aperture_cut_cell_mixer.h"

#include <functional>

#include "enums/interface_tag_definition.h"
#include "levelset/multi_phase_manager/two_phase_manager.h"

namespace {
constexpr unsigned int aperture_mixer_number_of_mixing_contributions = 9;
}

/**
 * @brief The default constructor of the ApertureCutCellMixer class. Calls the
 * default constructor of the base class.
 * @param halo_manager Instance to a HaloManager which provides MPI-related
 * methods.
 */
ApertureCutCellMixer::ApertureCutCellMixer(
    HaloManager &halo_manager, MaterialManager const &material_manager)
    : TwoPhaseCutCellMixer(halo_manager,
                           aperture_mixer_number_of_mixing_contributions,
                           material_manager) {
  // Empty Constructor, besides call of base class constructor
}

/**
 * @brief Determines all cell pairs which are mixed. One cell pair that is mixed
 * is denoted as mixing contribution.
 * @param node The node for which the mixing contributions are determined.
 * @param material The material which allows to identify the phase for which the
 * mixing contributions are determined.
 * @param mixing_contributions Indirect return which contains information about
 * the single mixing contributions. A single mixing contribution contains
 * information about the target cell indices (i_target, j_target and k_target -
 * saved as unsigned int), the mixing fraction beta (saved as double) and a
 * factor necessary to calculate the mixing fluxes.
 */
void ApertureCutCellMixer::CalculateMixingContributionsImplementation(
    Node const &node, MaterialName const material,
    std::vector<std::pair<std::vector<std::array<unsigned int, 6>>,
                          std::vector<std::array<double, 2>>>>
        &mixing_contributions) const {

  std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
  double const(&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock().GetReinitializedBuffer(
          InterfaceDescription::Levelset);
  double const(&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock().GetReinitializedBuffer(
          InterfaceDescription::VolumeFraction);
  std::int8_t const material_sign =
      MaterialSignCapsule::SignOfMaterial(material);

  double const reference_volume_fraction = (material_sign > 0) ? 0.0 : 1.0;
  double const material_sign_double = double(material_sign);

  // quantities to specify a mixing operation
  double volume_fraction_self = 0.0;
  double volume_fraction_target = 0.0;
  double mixing_fraction = 0.0;
  unsigned int i_target = 0;
  unsigned int j_target = 0;
  unsigned int k_target = 0;
  double beta_sum = 0.0; // temporary double to normalize mixing fractions
  bool mixing_for_cell_active = false;

  // vectors to describe the mixing operations for one cell
  std::vector<std::array<unsigned int, 6>> mixing_indices;
  std::vector<std::array<double, 2>> mixing_flux_factors;

  // lambda function to create a single mixing contribution between two cells
  // JW to HoNi This lambda function might be a performance bottle-neck
  std::function<void(unsigned int, unsigned int, unsigned int)>
      create_mixing_contribution =
          [&](unsigned int i, unsigned int j, unsigned int k) {
            double const denominator =
                volume_fraction_self * mixing_fraction + volume_fraction_target;
            if (mixing_fraction != 0.0 && denominator != 0.0) {
              mixing_for_cell_active = true;
              beta_sum += mixing_fraction;
              mixing_indices.push_back({i, j, k, i_target, j_target, k_target});
              mixing_flux_factors.push_back(
                  {mixing_fraction, volume_fraction_target});
            }
          };

  for (unsigned int i = FICMOX; i <= LICPOX; ++i) {
    for (unsigned int j = FICMOY; j <= LICPOY; ++j) {
      for (unsigned int k = FICMOZ; k <= LICPOZ; ++k) {
        // reset variables to prepare the calculations for the next cell
        mixing_indices.clear();
        mixing_flux_factors.clear();
        beta_sum = 0.0;
        mixing_for_cell_active = false;

        volume_fraction_self = reference_volume_fraction +
                               material_sign_double * volume_fraction[i][j][k];
        if (std::abs(interface_tags[i][j][k]) <= ITTI(IT::CutCellNeighbor) &&
            (volume_fraction_self < CC::MITH() ||
             levelset[i][j][k] * material_sign < 0)) {

          if (std::abs(interface_tags[i][j][k]) <= ITTI(IT::NewCutCell) &&
              volume_fraction_self != 0.0) {

            std::array<double, 6> const cell_face_apertures =
                geometry_calculator_.ComputeCellFaceAperture(levelset, i, j, k,
                                                             material_sign);

            // x-1
            i_target = i - 1;
            j_target = j;
            k_target = k;
            volume_fraction_target =
                reference_volume_fraction +
                material_sign_double *
                    volume_fraction[i_target][j_target][k_target];
            if (volume_fraction_target > volume_fraction_self) {
              mixing_fraction = std::abs(cell_face_apertures[0]);
              create_mixing_contribution(i, j, k);
            }

            // x+1
            i_target = i + 1;
            j_target = j;
            k_target = k;
            volume_fraction_target =
                reference_volume_fraction +
                material_sign_double *
                    volume_fraction[i_target][j_target][k_target];
            if (volume_fraction_target > volume_fraction_self) {
              mixing_fraction = std::abs(cell_face_apertures[1]);
              create_mixing_contribution(i, j, k);
            }

            if constexpr (CC::DIM() != Dimension::One) {
              // y-1
              i_target = i;
              j_target = j - 1;
              k_target = k;
              volume_fraction_target =
                  reference_volume_fraction +
                  material_sign_double *
                      volume_fraction[i_target][j_target][k_target];
              if (volume_fraction_target > volume_fraction_self) {
                mixing_fraction = std::abs(cell_face_apertures[2]);
                create_mixing_contribution(i, j, k);
              }

              // y+1
              i_target = i;
              j_target = j + 1;
              k_target = k;
              volume_fraction_target =
                  reference_volume_fraction +
                  material_sign_double *
                      volume_fraction[i_target][j_target][k_target];
              if (volume_fraction_target > volume_fraction_self) {
                mixing_fraction = std::abs(cell_face_apertures[3]);
                create_mixing_contribution(i, j, k);
              }
            }

            if constexpr (CC::DIM() == Dimension::Three) {
              // z-1
              i_target = i;
              j_target = j;
              k_target = k - 1;
              volume_fraction_target =
                  reference_volume_fraction +
                  material_sign_double *
                      volume_fraction[i_target][j_target][k_target];
              if (volume_fraction_target > volume_fraction_self) {
                mixing_fraction = std::abs(cell_face_apertures[4]);
                create_mixing_contribution(i, j, k);
              }

              // z+1
              i_target = i;
              j_target = j;
              k_target = k + 1;
              volume_fraction_target =
                  reference_volume_fraction +
                  material_sign_double *
                      volume_fraction[i_target][j_target][k_target];
              if (volume_fraction_target > volume_fraction_self) {
                mixing_fraction = std::abs(cell_face_apertures[5]);
                create_mixing_contribution(i, j, k);
              }
            }

            // normalization of the mixing fraction
            if (mixing_for_cell_active) {
              double const one_beta_sum = 1.0 / beta_sum;
              for (unsigned int n = 0; n < mixing_flux_factors.size(); ++n) {
                mixing_flux_factors[n][0] *= one_beta_sum;
              }
            }

          } else {

            mixing_for_cell_active = true;
            // calculate inward pointing normal vector
            std::array<double, 3> const normal =
                GetNormal(levelset, i, j, k, material_sign);

            // target cell in x direction
            i_target =
                i +
                ((normal[0] > 0.0) -
                 (normal[0] <
                  0.0)); // Bool arithmetics on purpose. +1 if normal direction
                         // positive, -1 if normal direction negative
            j_target = j;
            k_target = k;
            volume_fraction_target =
                reference_volume_fraction +
                material_sign_double *
                    volume_fraction[i_target][j_target][k_target];
            if (volume_fraction_target > volume_fraction_self) {
              mixing_fraction = normal[0] * normal[0];
              create_mixing_contribution(i, j, k);
            }

            if constexpr (CC::DIM() != Dimension::One) {
              // target cell in y direction
              i_target = i;
              j_target =
                  j + ((normal[1] > 0.0) -
                       (normal[1] < 0.0)); // Bool arithmetics on purpose. +1 if
                                           // normal direction positive, -1 if
                                           // normal direction negative
              k_target = k;
              volume_fraction_target =
                  reference_volume_fraction +
                  material_sign_double *
                      volume_fraction[i_target][j_target][k_target];
              if (volume_fraction_target > volume_fraction_self) {
                mixing_fraction = normal[1] * normal[1];
                create_mixing_contribution(i, j, k);
              }
            }

            if constexpr (CC::DIM() == Dimension::Three) {
              // target cell in z direction
              i_target = i;
              j_target = j;
              k_target =
                  k + ((normal[2] > 0.0) -
                       (normal[2] < 0.0)); // Bool arithmetics on purpose. +1 if
                                           // normal direction positive, -1 if
                                           // normal direction negative
              volume_fraction_target =
                  reference_volume_fraction +
                  material_sign_double *
                      volume_fraction[i_target][j_target][k_target];
              if (volume_fraction_target > volume_fraction_self) {
                mixing_fraction = normal[2] * normal[2];
                create_mixing_contribution(i, j, k);
              }
            }
          }

          if (mixing_for_cell_active) {
            // calculation of a constant to calculate the mixing fluxes
            for (unsigned int n = 0; n < mixing_flux_factors.size(); ++n) {
              // Note, that at this point in mixing_flux_factors[n][1] the
              // volume fraction is saved temporarily. This is necessary, since
              // the calculation of the mixing flux factor is only possible with
              // normalized mixing_fraction (saved in mixing_flux_factors[n][0])
              mixing_flux_factors[n][1] =
                  mixing_flux_factors[n][0] /
                  (volume_fraction_self * mixing_flux_factors[n][0] +
                   mixing_flux_factors[n][1]);
            }

            // add the mixing operations for cell i j k to the data container
            // which saves all mixing operations for one block
            mixing_contributions.push_back(
                std::make_pair(mixing_indices, mixing_flux_factors));
          }

        } // cells which are mixed
      }   // k
    }     // j
  }       // i
}
