//===------------------------ viscous_fluxes.cpp --------------------------===//
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
#include "viscous_fluxes.h"

#include <numeric>

#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "stencils/stencil_utilities.h"
#include "utilities/buffer_operations_stencils.h"
#include "utilities/index_transformations.h"

/**
 * @brief The constructor for the ViscousFluxes class.
 * @param material_manager The material manager which contains information about
 * the viscosities.
 */
ViscousFluxes::ViscousFluxes(MaterialManager const &material_manager)
    : material_manager_(material_manager) {
  // Empty constructor besides initializer list
}

/**
 * @brief Computes cell face fluxes due to viscosity.
 * @param mat_block A pair containing a block and its material.
 * @param dissipative_flux_x, dissipative_flux_y, dissipative_flux_z Reference
 * to the face fluxes (indirect return parameter).
 * @param cell_size The cell size.
 */
void ViscousFluxes::ComputeFluxes(
    std::pair<MaterialName const, Block> const &mat_block,
    double (&dissipative_flux_x)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                                [CC::ICZ() + 1],
    double (&dissipative_flux_y)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                                [CC::ICZ() + 1],
    double (&dissipative_flux_z)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1]
                                [CC::ICZ() + 1],
    double cell_size) const {

  // Define concretizations of the different stencil computations
  using DerivativeStencilCenter = DerivativeStencilSetup::Concretize<
      viscous_fluxes_derivative_stencil_cell_center>::type;
  using DerivativeStencilFace = DerivativeStencilSetup::Concretize<
      viscous_fluxes_derivative_stencil_cell_face>::type;
  using ReconstructionStencil = ReconstructionStencilSetup::Concretize<
      viscous_fluxes_reconstruction_stencil>::type;

  // y and z velocity buffers may not be available
  // as workaround use the x velocity buffer in these cases
  // this is legal since the respective gradients are not computed/used anyway
  double const(&u)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      mat_block.second.GetPrimeStateBuffer(PrimeState::VelocityX);
  double const(&v)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      CC::DIM() != Dimension::One
          ? mat_block.second.GetPrimeStateBuffer(PrimeState::VelocityY)
          : u;
  double const(&w)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      CC::DIM() == Dimension::Three
          ? mat_block.second.GetPrimeStateBuffer(PrimeState::VelocityZ)
          : u;

  // Get the shear viscosity from the material (all computations below ensure
  // that buffer is only used when model is active, otherwirse the buffer does
  // not exist)
  double const(&shear_viscosity)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      mat_block.second.GetParameterBuffer(Parameter::ShearViscosity);

  /**
   * Description for the positions of the Array:
   * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]
   * [DTI(CC::DIM())][DTI(CC::DIM())] Field index x  Field index y  Field index
   * z  Cell face x/y/z   Velocity gradient: du_i / dx_j
   */
  double velocity_gradient_at_cell_faces[CC::ICX() + 1][CC::ICY() + 1]
                                        [CC::ICZ() + 1][DTI(CC::DIM())]
                                        [DTI(CC::DIM())][DTI(CC::DIM())];

  /**
   * Description for the positions of the Array:
   * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]
   * [DTI(CC::DIM())] Field index x  Field index y  Field index z  Cell face
   * x/y/z   Velocity in x/y/z direction
   */
  double velocity_at_cell_faces[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                               [DTI(CC::DIM())][DTI(CC::DIM())];

  /**
   * Description for the positions of the Array:
   * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]
   * Field index x  Field index y  Field index z   Cell face x/y/z
   */
  double shear_viscosity_at_cell_faces[CC::ICX() + 1][CC::ICY() + 1]
                                      [CC::ICZ() + 1][DTI(CC::DIM())];

  /**
   * Description for the positions of the Array:
   * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]
   * [DTI(CC::DIM())][DTI(CC::DIM())] Field index x  Field index y  Field index
   * z  tau_ij at face i
   */
  double tau_flux[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI(CC::DIM())]
                 [DTI(CC::DIM())];

  for (unsigned int i = 0; i < CC::ICX() + 1; ++i) {
    for (unsigned int j = 0; j < CC::ICY() + 1; ++j) {
      for (unsigned int k = 0; k < CC::ICZ() + 1; ++k) {
        for (unsigned int r = 0; r < DTI(CC::DIM()); ++r) {
          for (unsigned int s = 0; s < DTI(CC::DIM()); ++s) {
            for (unsigned int t = 0; t < DTI(CC::DIM()); ++t) {
              velocity_gradient_at_cell_faces[i][j][k][r][s][t] = 0.0;
            }
            tau_flux[i][j][k][r][s] = 0.0;
            velocity_at_cell_faces[i][j][k][r][s] = 0.0;
          }

          if constexpr (CC::ShearViscosityModelActive()) {
            shear_viscosity_at_cell_faces[i][j][k][r] = 0.0;
          }
        }
      }
    }
  }
  // Compute tHe contributions for the dissipative fluxes
  BO::Stencils::ComputeVectorGradientAtCellFaces<
      DerivativeStencilCenter, DerivativeStencilFace, ReconstructionStencil>(
      u, v, w, cell_size, velocity_gradient_at_cell_faces);
  BO::Stencils::ComputeVectorAtCellFaces<ReconstructionStencil>(
      u, v, w, cell_size, velocity_at_cell_faces);
  // Depending if viscosity models are active, shear viscosity must be
  // reconstructed at cell face as well
  if constexpr (CC::ShearViscosityModelActive()) {
    BO::Stencils::ComputeScalarAtCellFaces<ReconstructionStencil>(
        shear_viscosity, cell_size, shear_viscosity_at_cell_faces);
    ComputeTauFluxes(
        velocity_gradient_at_cell_faces, shear_viscosity_at_cell_faces,
        material_manager_.GetMaterial(mat_block.first).GetBulkViscosity(),
        tau_flux);
  } else {
    ComputeTauFluxes(velocity_gradient_at_cell_faces,
                     material_manager_.GetMaterial(mat_block.first)
                         .GetShearAndBulkViscosity(),
                     tau_flux);
  }

  // Calculate the complete dissipative fluxes in each direction
  for (unsigned int i = 0; i < CC::ICX() + 1; ++i) {
    for (unsigned int j = 0; j < CC::ICY() + 1; ++j) {
      for (unsigned int k = 0; k < CC::ICZ() + 1; ++k) {

        for (unsigned int r = 0; r < DTI(CC::DIM()); ++r) {
          dissipative_flux_x[ETI(MF::AME()[r])][i][j][k] -=
              tau_flux[i][j][k][0][r];
          if constexpr (CC::DIM() != Dimension::One)
            dissipative_flux_y[ETI(MF::AME()[r])][i][j][k] -=
                tau_flux[i][j][k][1][r];
          if constexpr (CC::DIM() == Dimension::Three)
            dissipative_flux_z[ETI(MF::AME()[r])][i][j][k] -=
                tau_flux[i][j][k][2][r];
        }

        if constexpr (MF::IsEquationActive(Equation::Energy)) {
          double const energy_flux_x = std::inner_product(
              std::cbegin(tau_flux[i][j][k][0]),
              std::cend(tau_flux[i][j][k][0]),
              std::cbegin(velocity_at_cell_faces[i][j][k][0]), 0.0);
          dissipative_flux_x[ETI(Equation::Energy)][i][j][k] -= energy_flux_x;

          if constexpr (CC::DIM() != Dimension::One) {
            double const energy_flux_y = std::inner_product(
                std::cbegin(tau_flux[i][j][k][1]),
                std::cend(tau_flux[i][j][k][1]),
                std::cbegin(velocity_at_cell_faces[i][j][k][1]), 0.0);
            dissipative_flux_y[ETI(Equation::Energy)][i][j][k] -= energy_flux_y;
          }
          if constexpr (CC::DIM() == Dimension::Three) {
            double const energy_flux_z = std::inner_product(
                std::cbegin(tau_flux[i][j][k][2]),
                std::cend(tau_flux[i][j][k][2]),
                std::cbegin(velocity_at_cell_faces[i][j][k][2]), 0.0);
            dissipative_flux_z[ETI(Equation::Energy)][i][j][k] -= energy_flux_z;
          }
        }
      }
    }
  }
}

/**
 * @brief Computes tau, the viscous part of the stress tensor.
 * @param velocity_gradient_at_cell_faces The velocity gradient field at cell
 * faces.
 * @param viscosity A vector containing the shear and bulk viscosity.
 * @param tau The viscous part of the stress tensor.
 */
void ViscousFluxes::ComputeTauFluxes(
    double const (
        &velocity_gradient_at_cell_faces)[CC::ICX() + 1][CC::ICY() + 1]
                                         [CC::ICZ() + 1][DTI(CC::DIM())]
                                         [DTI(CC::DIM())][DTI(CC::DIM())],
    std::vector<double> const viscosity,
    double (&tau)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI(CC::DIM())]
                 [DTI(CC::DIM())]) const {

  double const mu_1 = viscosity[0];
  double const mu_2 = viscosity[1] - 2.0 * viscosity[0] / 3.0;

  for (unsigned int i = 0; i < CC::ICX() + 1; ++i) {
    for (unsigned int j = 0; j < CC::ICY() + 1; ++j) {
      for (unsigned int k = 0; k < CC::ICZ() + 1; ++k) {
        for (unsigned int r = 0; r < DTI(CC::DIM()); ++r) {
          double volumetric_part = 0.0;
          for (unsigned int s = 0; s < DTI(CC::DIM()); ++s) {
            tau[i][j][k][r][s] =
                mu_1 * (velocity_gradient_at_cell_faces[i][j][k][r][r][s] +
                        velocity_gradient_at_cell_faces[i][j][k][r][s][r]);
            volumetric_part +=
                velocity_gradient_at_cell_faces[i][j][k][r][s][s];
          }
          tau[i][j][k][r][r] += volumetric_part * mu_2;
        }
      }
    }
  }
}

/**
 * @brief Computes tau, the viscous part of the stress tensor, for a field
 * dependent viscosity.
 * @param velocity_gradient_at_cell_faces The velocity gradient field at cell
 * faces.
 * @param shear_viscosity The shear_viscosity field at the cell faces.
 * @param bulk_viscosity The constant bulk viscosity.
 * @param tau The viscous part of the stress tensor.
 */
void ViscousFluxes::ComputeTauFluxes(
    double const (
        &velocity_gradient_at_cell_faces)[CC::ICX() + 1][CC::ICY() + 1]
                                         [CC::ICZ() + 1][DTI(CC::DIM())]
                                         [DTI(CC::DIM())][DTI(CC::DIM())],
    double const (
        &shear_viscosity_at_cell_faces)[CC::ICX() + 1][CC::ICY() + 1]
                                       [CC::ICZ() + 1][DTI(CC::DIM())],
    double const bulk_viscosity,
    double (&tau)[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI(CC::DIM())]
                 [DTI(CC::DIM())]) const {

  for (unsigned int i = 0; i < CC::ICX() + 1; ++i) {
    for (unsigned int j = 0; j < CC::ICY() + 1; ++j) {
      for (unsigned int k = 0; k < CC::ICZ() + 1; ++k) {
        for (unsigned int r = 0; r < DTI(CC::DIM()); ++r) {
          // Compute the two viscosity coefficients
          double const mu_1 = shear_viscosity_at_cell_faces[i][j][k][r];
          double const mu_2 = bulk_viscosity - 2.0 * mu_1 / 3.0;

          double volumetric_part = 0.0;
          for (unsigned int s = 0; s < DTI(CC::DIM()); ++s) {
            tau[i][j][k][r][s] =
                mu_1 * (velocity_gradient_at_cell_faces[i][j][k][r][r][s] +
                        velocity_gradient_at_cell_faces[i][j][k][r][s][r]);
            volumetric_part +=
                velocity_gradient_at_cell_faces[i][j][k][r][s][s];
          }
          tau[i][j][k][r][r] += volumetric_part * mu_2;
        }
      }
    }
  }
}
