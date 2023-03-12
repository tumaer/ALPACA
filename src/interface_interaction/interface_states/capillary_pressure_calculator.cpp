//===-------------- capillary_pressure_calculator.cpp ---------------------===//
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
#include "capillary_pressure_calculator.h"
#include "enums/interface_tag_definition.h"
#include "levelset/geometry/geometry_calculator.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "stencils/stencil_utilities.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief      The default constructor of the class.
 *
 * @param[in]  material_manager  The material manager containing information
 * about the materials used in the simulation. Necessary to get the surface
 * tension coefficient.
 */
CapillaryPressureCalculator::CapillaryPressureCalculator(
    MaterialManager const &material_manager)
    : surface_tension_coefficient_(
          CC::CapillaryForcesActive()
              ? material_manager
                    .GetMaterialPairing(MaterialSignCapsule::NegativeMaterial(),
                                        MaterialSignCapsule::PositiveMaterial())
                    .GetSurfaceTensionCoefficient()
              : 0.0) {
  // Empty besides initializer list.
}

/**
 * @brief Calculates the capillary pressure in cut cells. Stores the result in a
 * buffer.
 *
 * @param node                 The node for which the capillary pressure is
 * calculated.
 * @param pressure_difference  The buffer in which the capillary pressure is
 * saved. Indirect return parameter.
 */
void CapillaryPressureCalculator::ComputePressureDifference(
    Node const &node,
    double (&pressure_difference)[CC::TCX()][CC::TCY()][CC::TCZ()]) const {

  std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
  /**
   * Use the pressure_difference buffer in the next function call to store the
   * curvature. Then, in the following step multiplication with the
   * surface_tension_coefficient_ gives the pressure difference due to capillary
   * forces.
   */
  ComputeInterfaceCurvature(node, pressure_difference);

  for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        if (std::abs(interface_tags[i][j][k]) <= ITTI(IT::NewCutCell)) {
          pressure_difference[i][j][k] *= surface_tension_coefficient_;
        }
      }
    }
  }
}
