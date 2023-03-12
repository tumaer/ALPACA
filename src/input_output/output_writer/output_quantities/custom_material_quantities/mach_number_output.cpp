//===-------------------- mach_number_output.cpp --------------------------===//
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
#include "input_output/output_writer/output_quantities/custom_material_quantities/mach_number_output.h"

#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "stencils/stencil_utilities.h"
#include "utilities/mathematical_functions.h"
#include "utilities/vector_utilities.h"
#include <algorithm> //lower_bound, sort

/**
 * @brief constructor to create the material Mach number output.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @param material_manager Instance to access all material data.
 * @param quantity_name Name of the quantity that is displayed in the ParaView
 * cell data list.
 * @param output_flags Flags of the output type that is written (0: standard, 1:
 * interface, 2:debug).
 *
 * @note {row, colmun} = {1,1} marks that the quantity is a scalar.
 */
MachNumberOutput::MachNumberOutput(UnitHandler const &unit_handler,
                                   MaterialManager const &material_manager,
                                   std::string const &quantity_name,
                                   std::array<bool, 3> const output_flags)
    : OutputQuantity(unit_handler, material_manager, quantity_name,
                     output_flags, {1, 1}) {
  /** Empty besides initializer list */
}

/**
 * @brief see base class definition.
 */
void MachNumberOutput::DoComputeCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter) const {

  if (node.HasLevelset()) {
    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

    PrimeStates const &positive_prime_states =
        node.GetPhaseByMaterial(MaterialSignCapsule::PositiveMaterial())
            .GetPrimeStateBuffer();
    PrimeStates const &negative_prime_states =
        node.GetPhaseByMaterial(MaterialSignCapsule::NegativeMaterial())
            .GetPrimeStateBuffer();

    EquationOfState const &positive_material_eos =
        material_manager_.GetMaterial(MaterialSignCapsule::PositiveMaterial())
            .GetEquationOfState();
    EquationOfState const &negative_material_eos =
        material_manager_.GetMaterial(MaterialSignCapsule::NegativeMaterial())
            .GetEquationOfState();

    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          if (interface_tags[i][j][k] > 0) {
            cell_data[cell_data_counter++] =
                VU::L2Norm(
                    {positive_prime_states[PrimeState::VelocityX][i][j][k],
                     CC::DIM() != Dimension::One
                         ? positive_prime_states[PrimeState::VelocityY][i][j][k]
                         : 0.0,
                     CC::DIM() == Dimension::Three
                         ? positive_prime_states[PrimeState::VelocityZ][i][j][k]
                         : 0.0}) /
                positive_material_eos.SpeedOfSound(
                    positive_prime_states[PrimeState::Density][i][j][k],
                    positive_prime_states[PrimeState::Pressure][i][j][k]);
          } else {
            cell_data[cell_data_counter++] =
                VU::L2Norm(
                    {negative_prime_states[PrimeState::VelocityX][i][j][k],
                     CC::DIM() != Dimension::One
                         ? negative_prime_states[PrimeState::VelocityY][i][j][k]
                         : 0.0,
                     CC::DIM() == Dimension::Three
                         ? negative_prime_states[PrimeState::VelocityZ][i][j][k]
                         : 0.0}) /
                negative_material_eos.SpeedOfSound(
                    negative_prime_states[PrimeState::Density][i][j][k],
                    negative_prime_states[PrimeState::Pressure][i][j][k]);
          }
        }
      }
    }
  } else {
    // No interface node -> interface tags/material is the same everywhere
    MaterialName const material = node.GetSinglePhaseMaterial();
    PrimeStates const &prime_states =
        node.GetPhaseByMaterial(material).GetPrimeStateBuffer();
    EquationOfState const &material_eos =
        material_manager_.GetMaterial(material).GetEquationOfState();

    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          cell_data[cell_data_counter++] =
              VU::L2Norm({prime_states[PrimeState::VelocityX][i][j][k],
                          CC::DIM() != Dimension::One
                              ? prime_states[PrimeState::VelocityY][i][j][k]
                              : 0.0,
                          CC::DIM() == Dimension::Three
                              ? prime_states[PrimeState::VelocityZ][i][j][k]
                              : 0.0}) /
              material_eos.SpeedOfSound(
                  prime_states[PrimeState::Density][i][j][k],
                  prime_states[PrimeState::Pressure][i][j][k]);
        }
      }
    }
  }
}

/**
 * @brief see base class definition.
 *
 * @note Attention: In case prime state, parameter  variables are used, pay
 * attention that they only exist on leave nodes. In case a division is made on
 * non-leave nodes a floating point exception is caused. Therefore, only use the
 * debug output if it is ensured that this cannot happen. Conservatives can be
 * used since they are present on all nodes.
 */
void MachNumberOutput::DoComputeDebugCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter,
    MaterialName const material) const {

  // Compute the real value only for nodes that contain the material
  if (node.ContainsMaterial(material)) {
    for (unsigned int k = 0; k < CC::TCZ(); ++k) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int i = 0; i < CC::TCX(); ++i) {
          cell_data[cell_data_counter++] = 1.0;
        }
      }
    }
  } else {
    // otherwise use default value
    for (unsigned int k = 0; k < CC::TCZ(); ++k) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int i = 0; i < CC::TCX(); ++i) {
          cell_data[cell_data_counter++] = -1.0;
        }
      }
    }
  }
}
