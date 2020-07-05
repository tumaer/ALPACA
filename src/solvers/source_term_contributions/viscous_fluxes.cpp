/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#include "viscous_fluxes.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "utilities/index_transformations.h"
#include "stencils/stencil_utilities.h"
#include "utilities/buffer_operations_stencils.h"

/**
 * @brief The constructor for the ViscousFluxes class.
 * @param material_manager The material manager which contains information about the viscosities.
 */
ViscousFluxes::ViscousFluxes(MaterialManager const& material_manager) :
   material_manager_(material_manager)
{
   // Empty constructor besides initializer list
}

/**
 * @brief Computes cell face fluxes due to viscosity.
 * @param mat_block A pair containing a block and its material.
 * @param dissipative_flux_x, dissipative_flux_y, dissipative_flux_z Reference to the face fluxes (indirect return parameter).
 * @param cell_size The cell size.
 */
void ViscousFluxes::ComputeFluxes( std::pair<MaterialName const, Block> const& mat_block
                                 , double (&dissipative_flux_x)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                                 , double (&dissipative_flux_y)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                                 , double (&dissipative_flux_z)[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1]
                                 , double cell_size ) const {

   // Define concretizations of the different stencil computations
   using DerivativeStencilCenter = DerivativeStencilSetup::Concretize<viscous_fluxes_derivative_stencil_cell_center>::type;
   using DerivativeStencilFace   = DerivativeStencilSetup::Concretize<viscous_fluxes_derivative_stencil_cell_face>::type;
   using ReconstructionStencil   = ReconstructionStencilSetup::Concretize<viscous_fluxes_reconstruction_stencil>::type;

   // y and z velocity buffers may not be available
   // as workaround use the x velocity buffer in these cases
   // this is legal since the respective gradients are not computed/used anyway
   double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()] = mat_block.second.GetPrimeStateBuffer(PrimeState::VelocityX);
   double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()] = CC::DIM() != Dimension::One   ? mat_block.second.GetPrimeStateBuffer(PrimeState::VelocityY) : u;
   double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()] = CC::DIM() == Dimension::Three ? mat_block.second.GetPrimeStateBuffer(PrimeState::VelocityZ) : u;

   // Get the shear viscosity from the material (all computations below ensure that buffer is only used when model is active, otherwirse the buffer does not exist)
   double const (&shear_viscosity)[CC::TCX()][CC::TCY()][CC::TCZ()] = mat_block.second.GetParameterBuffer(Parameter::ShearViscosity);

   /**
    * Description for the positions of the Array:
    * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]  [DTI(CC::DIM())][DTI(CC::DIM())]
    * Field index x  Field index y  Field index z  Cell face x/y/z   Velocity gradient: du_i / dx_j
    */
   double velocity_gradient_at_cell_faces[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())][DTI(CC::DIM())];

   /**
    * Description for the positions of the Array:
    * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]  [DTI(CC::DIM())]
    * Field index x  Field index y  Field index z  Cell face x/y/z   Velocity in x/y/z direction
    */
   double velocity_at_cell_faces[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())];

   /**
    * Description for the positions of the Array:
    * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]
    * Field index x  Field index y  Field index z   Cell face x/y/z
    */
   double shear_viscosity_at_cell_faces[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())];

   /**
    * Description for the positions of the Array:
    * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())][DTI(CC::DIM())]
    * Field index x  Field index y  Field index z  tau_ij at face i
    */
   double tau_flux[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())];

   for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
      for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
         for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {
                  for( unsigned int t = 0; t < DTI(CC::DIM()); ++t ) {
                     velocity_gradient_at_cell_faces[i][j][k][r][s][t] = 0.0;
                  }
                  tau_flux[i][j][k][r][s] = 0.0;
                  velocity_at_cell_faces[i][j][k][r][s] = 0.0;
               }

               if constexpr( CC::ShearViscosityModelActive() ) {
                  shear_viscosity_at_cell_faces[i][j][k][r] = 0.0;
               }
            }
         }
      }
   }
   // Compute tHe contributions for the dissipative fluxes
   BO::Stencils::ComputeVectorGradientAtCellFaces<DerivativeStencilCenter,DerivativeStencilFace,ReconstructionStencil>(u,v,w,cell_size,velocity_gradient_at_cell_faces);
   BO::Stencils::ComputeVectorAtCellFaces<ReconstructionStencil>(u,v,w,cell_size,velocity_at_cell_faces);
   // Depending if viscosity models are active, shear viscosity must be reconstructed at cell face as well
   if constexpr( CC::ShearViscosityModelActive() ) {
      BO::Stencils::ComputeScalarAtCellFaces<ReconstructionStencil>(shear_viscosity, cell_size, shear_viscosity_at_cell_faces);
      ComputeTauFluxes(velocity_gradient_at_cell_faces, shear_viscosity_at_cell_faces, material_manager_.GetMaterial(mat_block.first).GetBulkViscosity(), tau_flux);
   }
   else {
      ComputeTauFluxes(velocity_gradient_at_cell_faces, material_manager_.GetMaterial(mat_block.first).GetShearAndBulkViscosity(), tau_flux);
   }

   // Calculate the complete dissipative fluxes in each direction
   for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
      for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
         for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
            double energy_flux_x = 0.0;
            double energy_flux_y = 0.0;
            double energy_flux_z = 0.0;

            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               dissipative_flux_x[ETI(MF::AME()[r])][i][j][k] -= tau_flux[i][j][k][0][r];
               if constexpr( CC::DIM() != Dimension::One )   dissipative_flux_y[ETI(MF::AME()[r])][i][j][k] -= tau_flux[i][j][k][1][r];
               if constexpr( CC::DIM() == Dimension::Three ) dissipative_flux_z[ETI(MF::AME()[r])][i][j][k] -= tau_flux[i][j][k][2][r];

               energy_flux_x += tau_flux[i][j][k][0][r] * velocity_at_cell_faces[i][j][k][0][r];
               if constexpr( CC::DIM() != Dimension::One )   energy_flux_y += tau_flux[i][j][k][1][r] * velocity_at_cell_faces[i][j][k][1][r];
               if constexpr( CC::DIM() == Dimension::Three ) energy_flux_z += tau_flux[i][j][k][2][r] * velocity_at_cell_faces[i][j][k][2][r];
            }

            dissipative_flux_x[ETI(Equation::Energy)][i][j][k] -= energy_flux_x;
            if constexpr( CC::DIM() != Dimension::One )   dissipative_flux_y[ETI(Equation::Energy)][i][j][k] -= energy_flux_y;
            if constexpr( CC::DIM() == Dimension::Three ) dissipative_flux_z[ETI(Equation::Energy)][i][j][k] -= energy_flux_z;

         }
      }
   }

}

/**
 * @brief Computes tau, the viscous part of the stress tensor.
 * @param velocity_gradient_at_cell_faces The velocity gradient field at cell faces.
 * @param viscosity A vector containing the shear and bulk viscosity.
 * @param tau The viscous part of the stress tensor.
 */
void ViscousFluxes::ComputeTauFluxes( double const (&velocity_gradient_at_cell_faces)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())][DTI(CC::DIM())]
                                     ,std::vector<double> const viscosity
                                     ,double (&tau)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())]) const {

   double const mu_1 = viscosity[0];
   double const mu_2 = viscosity[1] - 2.0 * viscosity[0] / 3.0;

   for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
      for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
         for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               double volumetric_part = 0.0;
               for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {
                  tau[i][j][k][r][s] = mu_1 * (velocity_gradient_at_cell_faces[i][j][k][r][r][s] + velocity_gradient_at_cell_faces[i][j][k][r][s][r]);
                  volumetric_part += velocity_gradient_at_cell_faces[i][j][k][r][s][s];
               }
               tau[i][j][k][r][r] += volumetric_part * mu_2;
            }
         }
      }
   }
}

/**
 * @brief Computes tau, the viscous part of the stress tensor, for a field dependent viscosity
 * @param velocity_gradient_at_cell_faces The velocity gradient field at cell faces.
 * @param shear_viscosity The shear_viscosity field at the cell faces
 * @param bulk_viscoisty The constant bulk viscosity
 * @param tau The viscous part of the stress tensor.
 */
void ViscousFluxes::ComputeTauFluxes( double const (&velocity_gradient_at_cell_faces)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())][DTI(CC::DIM())]
                                     ,double const (&shear_viscosity_at_cell_faces)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())]
                                     ,double const bulk_viscosity
                                     ,double (&tau)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][DTI(CC::DIM())][DTI(CC::DIM())]) const {

   for( unsigned int i = 0; i < CC::ICX() + 1; ++i ) {
      for( unsigned int j = 0; j < CC::ICY() + 1; ++j ) {
         for( unsigned int k = 0; k < CC::ICZ() + 1; ++k ) {
            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               // Compute the two viscosity coefficients
               double const mu_1 = shear_viscosity_at_cell_faces[i][j][k][r];
               double const mu_2 = bulk_viscosity - 2.0 * mu_1 / 3.0;

               double volumetric_part = 0.0;
               for( unsigned int s = 0; s < DTI(CC::DIM()); ++s ) {
                  tau[i][j][k][r][s] = mu_1 * (velocity_gradient_at_cell_faces[i][j][k][r][r][s] + velocity_gradient_at_cell_faces[i][j][k][r][s][r]);
                  volumetric_part += velocity_gradient_at_cell_faces[i][j][k][r][s][s];
               }
               tau[i][j][k][r][r] += volumetric_part * mu_2;
            }
         }
      }
   }
}