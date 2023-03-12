//===----------------- numerical_schlieren_output.cpp ---------------------===//
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
#include "input_output/output_writer/output_quantities/custom_material_quantities/numerical_schlieren_output.h"

#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "stencils/stencil_utilities.h"
#include "utilities/mathematical_functions.h"
#include "utilities/vector_utilities.h"
#include <algorithm> //lower_bound, sort

/**
 * @brief constructor to create a numerical schlieren output.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @param material_manager Instance to access all material data.
 * @param quantity_name Name of the quantity that is displayed in the ParaView
 * cell data list.
 * @param output_flags Flags of the output type that is written (0: standard, 1:
 * interface, 2:debug).
 *
 * @note {row, colmun} = {1,1} marks that the quantity is a scalar.
 */
NumericalSchlierenOutput::NumericalSchlierenOutput(
    UnitHandler const &unit_handler, MaterialManager const &material_manager,
    std::string const &quantity_name, std::array<bool, 3> const output_flags)
    : OutputQuantity(unit_handler, material_manager, quantity_name,
                     output_flags, {1, 1}) {
  /** Empty besides initializer list */
}

/**
 * @brief see base class definition.
 */
void NumericalSchlierenOutput::DoComputeCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter) const {
  // define derivative stencils for derivative computations
  constexpr DerivativeStencils derivative_stencil =
      DerivativeStencils::FourthOrderCentralDifference;
  using DerivativeStencil =
      DerivativeStencilSetup::Concretize<derivative_stencil>::type;

  // Obtain the factor used for dimensionalization unit: [kg/m^4]
  double const dimensionalization_factor = unit_handler_.DimensionalizeValue(
      1.0, {UnitType::Density}, {UnitType::Length});

  // get the cell size for the node
  double const cell_size = node.GetCellSize();

  if (node.HasLevelset()) {
    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
    double const(&positive_density)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetPhaseByMaterial(MaterialSignCapsule::PositiveMaterial())
            .GetPrimeStateBuffer(PrimeState::Density);
    double const(&negative_density)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetPhaseByMaterial(MaterialSignCapsule::NegativeMaterial())
            .GetPrimeStateBuffer(PrimeState::Density);
    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          if (interface_tags[i][j][k] > 0) {
            cell_data[cell_data_counter++] =
                VU::L2Norm(SU::GradientVector<DerivativeStencil>(
                    positive_density, i, j, k, cell_size)) *
                dimensionalization_factor;
          } else {
            cell_data[cell_data_counter++] =
                VU::L2Norm(SU::GradientVector<DerivativeStencil>(
                    negative_density, i, j, k, cell_size)) *
                dimensionalization_factor;
          }
        }
      }
    }
  } else {
    // No interface node -> interface tags/material is the same everywhere
    MaterialName const material = node.GetSinglePhaseMaterial();
    Block const &block = node.GetPhaseByMaterial(material);
    double const(&density)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::Density);
    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          cell_data[cell_data_counter++] =
              VU::L2Norm(SU::GradientVector<DerivativeStencil>(density, i, j, k,
                                                               cell_size)) *
              dimensionalization_factor;
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
void NumericalSchlierenOutput::DoComputeDebugCellData(
    Node const &node, std::vector<double> &cell_data,
    unsigned long long int &cell_data_counter,
    MaterialName const material) const {

  /** Now the actual assinging of to the hdf5 written data vector is done
   * depending on the given material */
  if (node.ContainsMaterial(material)) {
    for (unsigned int k = 0; k < CC::TCZ(); ++k) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int i = 0; i < CC::TCX(); ++i) {
          cell_data[cell_data_counter++] = 1.0;
        }
      }
    }
  } else {
    // otherwise use debug_default_value
    for (unsigned int k = 0; k < CC::TCZ(); ++k) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int i = 0; i < CC::TCX(); ++i) {
          cell_data[cell_data_counter++] = -1.0;
        }
      }
    }
  }
}
