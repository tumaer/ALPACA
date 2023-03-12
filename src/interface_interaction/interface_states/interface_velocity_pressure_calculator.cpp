//===--------- interface_velocity_pressure_calculator.cpp -----------------===//
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
#include "interface_velocity_pressure_calculator.h"
#include "enums/interface_tag_definition.h"
#include "levelset/geometry/geometry_calculator.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief      Default constructor for the InterfaceVelocityPressureCalculator
 * class. Initializes the members.
 *
 * @param[in]  material_manager  The material manager containing information
 * about the materials involved in the simulation.
 */
InterfaceVelocityPressureCalculator::InterfaceVelocityPressureCalculator(
    MaterialManager const &material_manager)
    : interface_riemann_solver_(material_manager) {
  // Empty besides initializer list
}

/**
 * @brief      Fills the interface velocity and pressure buffer of the levelset
 * block of a given node. This function may consider a pressure jump due to
 * capillary forces at the interface.
 *
 * @param      node                 The node for which the interface velocity
 * and pressure are calculated.
 * @param      pressure_difference  The pressure difference induced by capillary
 * forces.
 */
void InterfaceVelocityPressureCalculator::
    FillInterfaceVelocityAndPressureBuffer(
        Node &node,
        double (&pressure_difference)[CC::TCX()][CC::TCY()][CC::TCZ()]) const {

  double const(&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock().GetReinitializedBuffer(
          InterfaceDescription::Levelset);
  std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

  double(&interface_velocity)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock().GetInterfaceStateBuffer(
          InterfaceState::Velocity);
  double(&interface_pressure_positive)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock().GetInterfaceStateBuffer(
          InterfaceState::PressurePositive);
  double(&interface_pressure_negative)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      CC::CapillaryForcesActive()
          ? node.GetInterfaceBlock().GetInterfaceStateBuffer(
                InterfaceState::PressureNegative)
          : node.GetInterfaceBlock().GetInterfaceStateBuffer(
                InterfaceState::PressurePositive);

  // for describing positive material as right, and negative material als left,
  // see cited paper of Luo
  MaterialName const material_left = MaterialSignCapsule::NegativeMaterial();
  MaterialName const material_right = MaterialSignCapsule::PositiveMaterial();

  PrimeStates const &left_prime_states =
      node.GetPhaseByMaterial(material_left).GetPrimeStateBuffer();
  PrimeStates const &right_prime_states =
      node.GetPhaseByMaterial(material_right).GetPrimeStateBuffer();

  double velocity_normal_left = 0.0;
  double velocity_normal_right = 0.0;

  for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        if (std::abs(interface_tags[i][j][k]) <= ITTI(IT::ExtensionBand)) {

          std::array<double, 3> const normal =
              GetNormal(levelset_reinitialized, i, j, k);

          velocity_normal_left =
              left_prime_states[PrimeState::VelocityX][i][j][k] * normal[0];
          velocity_normal_left +=
              CC::DIM() != Dimension::One
                  ? left_prime_states[PrimeState::VelocityY][i][j][k] *
                        normal[1]
                  : 0.0;
          velocity_normal_left +=
              CC::DIM() == Dimension::Three
                  ? left_prime_states[PrimeState::VelocityZ][i][j][k] *
                        normal[2]
                  : 0.0;

          velocity_normal_right =
              right_prime_states[PrimeState::VelocityX][i][j][k] * normal[0];
          velocity_normal_right +=
              CC::DIM() != Dimension::One
                  ? right_prime_states[PrimeState::VelocityY][i][j][k] *
                        normal[1]
                  : 0.0;
          velocity_normal_right +=
              CC::DIM() == Dimension::Three
                  ? right_prime_states[PrimeState::VelocityZ][i][j][k] *
                        normal[2]
                  : 0.0;

          std::array<double, 3> const interface_states =
              interface_riemann_solver_.SolveInterfaceRiemannProblem(
                  left_prime_states[PrimeState::Density][i][j][k],
                  left_prime_states[PrimeState::Pressure][i][j][k],
                  velocity_normal_left, material_left,
                  right_prime_states[PrimeState::Density][i][j][k],
                  right_prime_states[PrimeState::Pressure][i][j][k],
                  velocity_normal_right, material_right,
                  pressure_difference[i][j][k]);

          interface_velocity[i][j][k] = interface_states[0];
          interface_pressure_positive[i][j][k] = interface_states[1];
          if constexpr (CC::CapillaryForcesActive()) {
            interface_pressure_negative[i][j][k] = interface_states[2];
          }
        } // if
      }   // k
    }     // j
  }       // i
}
