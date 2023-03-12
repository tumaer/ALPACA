//===------------------------- heat_fluxes.cpp ----------------------------===//
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
#include "heat_fluxes.h"
#include "stencils/stencil_utilities.h"
#include "utilities/buffer_operations_stencils.h"

/**
 * @brief The default constructor for the class.
 * @param material_manager The material manager.
 */
HeatFluxes::HeatFluxes(MaterialManager const &material_manager)
    : material_manager_(material_manager) {
  // Empty constructor besides initializer list
}

/**
 * @brief Computes the fluxes and adds them to buffers for fluxes in x-, y- and
 * z-direction.
 * @param mat_block The material block pair of the fluid for which the fluxes
 * are calculated.
 * @param heat_fluxes_x The heat fluxes in x-direction (indirect return
 * parameter).
 * @param heat_fluxes_y The heat fluxes in y-direction (indirect return
 * parameter).
 * @param heat_fluxes_z The heat fluxes in z-direction (indirect return
 * parameter).
 * @param cell_size The cell size of the block.
 */
void HeatFluxes::ComputeFluxes(
    std::pair<MaterialName const, Block> const &mat_block,
    double (&heat_fluxes_x)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                           [CC::ICZ() + 1],
    double (&heat_fluxes_y)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                           [CC::ICZ() + 1],
    double (&heat_fluxes_z)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                           [CC::ICZ() + 1],
    double const cell_size) const {

  // Definition of all stencils used during the flux calculation
  using DerivativeStencilCellFace = DerivativeStencilSetup::Concretize<
      heat_fluxes_derivative_stencil_cell_face>::type;
  using ReconstructionStencil = ReconstructionStencilSetup::Concretize<
      heat_fluxes_reconstruction_stencil>::type;

  // Fixed thermal conductivity (only used if no model active, otherwise also
  // non meaningful value)
  double const thermal_conductivity =
      material_manager_.GetMaterial(mat_block.first).GetThermalConductivity();

  // Get the temperature field buffer
  double const(&temperature)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      mat_block.second.GetPrimeStateBuffer(PrimeState::Temperature);

  // Get the thermal conductivity from the material (all computations below
  // ensure that buffer is only used when model is active, otherwise the buffer
  // does not exist)
  double const(&conductivity)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      mat_block.second.GetParameterBuffer(Parameter::ThermalConductivity);

  /**
   * Description for the positions of the Array:
   * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]
   * Field index x  Field index y  Field index z   Cell face x/y/z
   */
  double conductivity_at_cell_faces[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                                   [DTI(CC::DIM())];

  // Only fill the buffer if model is active
  if constexpr (CC::ThermalConductivityModelActive()) {
    for (unsigned int i = 0; i < CC::ICX() + 1; ++i) {
      for (unsigned int j = 0; j < CC::ICY() + 1; ++j) {
        for (unsigned int k = 0; k < CC::ICZ() + 1; ++k) {
          for (unsigned int r = 0; r < DTI(CC::DIM()); ++r) {
            conductivity_at_cell_faces[i][j][k][r] = 0.0;
          }
        }
      }
    }

    // Compute the conductivity properly
    BO::Stencils::ComputeScalarAtCellFaces<ReconstructionStencil>(
        conductivity, cell_size, conductivity_at_cell_faces);
  }

  // offset for flux-arrays
  constexpr int offset_x = CC::FICX() - 1;
  constexpr int offset_y = CC::DIM() != Dimension::One ? CC::FICY() - 1 : -1;
  constexpr int offset_z = CC::DIM() == Dimension::Three ? CC::FICZ() - 1 : -1;

  // Compute changes due to heat transfer - compute second order derivative of
  // temperature with fourth-order central-difference stencil
  for (unsigned int i = CC::FICX() - 1; i <= CC::LICX(); ++i) {
    for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        if constexpr (CC::ThermalConductivityModelActive()) {
          // take the thermal conductivity at the x cell face ( 0 index )
          heat_fluxes_x[ETI(Equation::Energy)][i - offset_x][j - offset_y]
                       [k - offset_z] +=
              -conductivity_at_cell_faces[i][j][k][0] *
              SU::Derivative<DerivativeStencilCellFace, Direction::X>(
                  temperature, i, j, k, cell_size);
        } else {
          heat_fluxes_x[ETI(Equation::Energy)][i - offset_x][j - offset_y]
                       [k - offset_z] +=
              -thermal_conductivity *
              SU::Derivative<DerivativeStencilCellFace, Direction::X>(
                  temperature, i, j, k, cell_size);
        }
      }
    }
  }

  if constexpr (CC::DIM() != Dimension::One) {
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY() - 1; j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
          if constexpr (CC::ThermalConductivityModelActive()) {
            // take the thermal conductivity at the y cell face ( 1 index )
            heat_fluxes_y[ETI(Equation::Energy)][i - offset_x][j - offset_y]
                         [k - offset_z] +=
                -conductivity_at_cell_faces[i][j][k][1] *
                SU::Derivative<DerivativeStencilCellFace, Direction::Y>(
                    temperature, i, j, k, cell_size);
          } else {
            heat_fluxes_y[ETI(Equation::Energy)][i - offset_x][j - offset_y]
                         [k - offset_z] +=
                -thermal_conductivity *
                SU::Derivative<DerivativeStencilCellFace, Direction::Y>(
                    temperature, i, j, k, cell_size);
          }
        }
      }
    }
  }

  if constexpr (CC::DIM() == Dimension::Three) {
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ() - 1; k <= CC::LICZ(); ++k) {
          if constexpr (CC::ThermalConductivityModelActive()) {
            // take the thermal conductivity at the z cell face ( 2 index )
            heat_fluxes_z[ETI(Equation::Energy)][i - offset_x][j - offset_y]
                         [k - offset_z] +=
                -conductivity_at_cell_faces[i][j][k][2] *
                SU::Derivative<DerivativeStencilCellFace, Direction::Z>(
                    temperature, i, j, k, cell_size);
          } else {
            heat_fluxes_z[ETI(Equation::Energy)][i - offset_x][j - offset_y]
                         [k - offset_z] +=
                -thermal_conductivity *
                SU::Derivative<DerivativeStencilCellFace, Direction::Z>(
                    temperature, i, j, k, cell_size);
          }
        }
      }
    }
  }
}
