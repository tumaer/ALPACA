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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#include "interface_stress_tensor_fluxes.h"
#include "utilities/mathematical_functions.h"
#include "utilities/index_transformations.h"
#include "levelset/geometry/geometry_calculator.h"
#include "enums/interface_tag_definition.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief      Default constructor of the class.
 *
 * @param[in]  material_positive  The material of the positive material.
 * @param[in]  mu_positive        The shear and bulk viscosity of the positive material.
 * @param[in]  material_negative  The material of the negative material.
 * @param[in]  mu_negative        The shear and bulk viscosity of the negative material.
 */
InterfaceStressTensorFluxes::InterfaceStressTensorFluxes( MaterialManager const& material_manager, MaterialName const material_positive, std::vector<double> const mu_positive,
                                                          MaterialName const material_negative, std::vector<double> const mu_negative ) : material_manager_( material_manager ),
                                                                                                                                          positive_material_properties_( material_positive, mu_positive ),
                                                                                                                                          negative_material_properties_( material_negative, mu_negative ) {
}

/**
 * @brief      Calculates the interface viscosities, i.e. the shear, bulk and volumetric contribution, by a volumetric averaging.
 *
 * @param[in]  volume_fraction  The volume fraction used for the averaging.
 *
 * @return     The shear, bulk and volumetric viscosities comprised in a vector.
 */
std::array<double, 3> InterfaceStressTensorFluxes::ComputeInterfaceViscosities( double const volume_fraction ) const {
   double const mu_shear_interface = positive_material_properties_.mu_shear_ * negative_material_properties_.mu_shear_ / ( volume_fraction * negative_material_properties_.mu_shear_ + ( 1 - volume_fraction ) * positive_material_properties_.mu_shear_ + epsilon_ );
   double const mu_bulk_interface  = positive_material_properties_.mu_bulk_ * negative_material_properties_.mu_bulk_ / ( volume_fraction * negative_material_properties_.mu_bulk_ + ( 1 - volume_fraction ) * positive_material_properties_.mu_bulk_ + epsilon_ );
   double const mu_lame_interface  = mu_bulk_interface - 2.0 * mu_shear_interface / 3.0;
   return { mu_shear_interface, mu_bulk_interface, mu_lame_interface };
}

/**
 * @brief      Computes the interface fluxes due to the stress tensor and adds them to the right-hand side buffer.
 *
 * @param      node                      The node for which the fluxes are calculated.
 * @param[in]  delta_aperture_field      The delta aperture field.
 * @param[in]  u_interface_normal_field  The field of the normal projected interface velocity.
 */
void InterfaceStressTensorFluxes::ComputeInterfaceFluxes( Node& node, double const ( &delta_aperture_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3], double const ( &u_interface_normal_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3] ) const {

   double interface_stress_positive[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )];
   double interface_stress_negative[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )];
   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
               for( unsigned int s = 0; s < DTI( CC::DIM() ); ++s ) {
                  interface_stress_positive[i][j][k][r][s] = 0.0;
                  interface_stress_negative[i][j][k][r][s] = 0.0;
               }// s
            }   // r
         }      // k
      }         // j
   }            // i

   if( CC::InviscidExchangeActive() ) {
      AddInviscidPartToInterfaceStressTensor( node, interface_stress_positive, interface_stress_negative );
   }

   if( CC::ViscosityIsActive() ) {
      AddViscousPartToInterfaceStressTensor( node, interface_stress_positive, interface_stress_negative );
   }

   AddFluxesToRightHandSide( node, delta_aperture_field, u_interface_normal_field, interface_stress_positive, interface_stress_negative );
}

/**
 * @brief      Adds the fluxes to the right-hand side buffer.
 *
 * @param      node                                    The node for which the fluxes are added.
 * @param[in]  delta_aperture_field                    The delta aperture field.
 * @param[in]  u_interface_normal_field                The normal projected interface velocity field.
 * @param[in]  interface_stress_tensor_positive_material  The interface stress tensor of the positive material.
 * @param[in]  interface_stress_tensor_negative_material  The interface stress tensor of the negative material.
 */
void InterfaceStressTensorFluxes::AddFluxesToRightHandSide( Node& node, double const ( &delta_aperture_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3], double const ( &u_interface_normal_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3], double const ( &interface_stress_tensor_positive_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )], double const ( &interface_stress_tensor_negative_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const {

   double const one_cell_size                                            = 1.0 / node.GetCellSize();
   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

   Conservatives& right_hand_side_positive_material = node.GetPhaseByMaterial( positive_material_properties_.material_ ).GetRightHandSideBuffer();
   Conservatives& right_hand_side_negative_material = node.GetPhaseByMaterial( negative_material_properties_.material_ ).GetRightHandSideBuffer();

   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            std::array<unsigned int, 3> const indices = { BIT::I2TX( i ), BIT::I2TY( j ), BIT::I2TZ( k ) };

            if( std::abs( interface_tags[indices[0]][indices[1]][indices[2]] ) <= ITTI( IT::NewCutCell ) ) {

               std::array<double, DTI( CC::DIM() )> momentum_fluxes_positive_material;
               std::array<double, DTI( CC::DIM() )> momentum_fluxes_negative_material;
               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  momentum_fluxes_positive_material[r] = 0.0;
                  momentum_fluxes_negative_material[r] = 0.0;
                  for( unsigned int s = 0; s < DTI( CC::DIM() ); ++s ) {
                     momentum_fluxes_positive_material[r] += interface_stress_tensor_positive_material[i][j][k][r][s] * delta_aperture_field[i][j][k][s];
                     momentum_fluxes_negative_material[r] += interface_stress_tensor_negative_material[i][j][k][r][s] * delta_aperture_field[i][j][k][s];
                  }
               }

               double enery_flux_positive_material = 0.0;
               double enery_flux_negative_material = 0.0;
               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  enery_flux_positive_material += momentum_fluxes_positive_material[r] * u_interface_normal_field[i][j][k][r];
                  enery_flux_negative_material += momentum_fluxes_negative_material[r] * u_interface_normal_field[i][j][k][r];
               }
               if constexpr( MF::IsEquationActive( Equation::Energy ) ) {
                  right_hand_side_positive_material[Equation::Energy][indices[0]][indices[1]][indices[2]] -= enery_flux_positive_material * one_cell_size;
                  right_hand_side_negative_material[Equation::Energy][indices[0]][indices[1]][indices[2]] += enery_flux_negative_material * one_cell_size;
               }
               right_hand_side_positive_material[Equation::MomentumX][indices[0]][indices[1]][indices[2]] -= momentum_fluxes_positive_material[0] * one_cell_size;
               right_hand_side_negative_material[Equation::MomentumX][indices[0]][indices[1]][indices[2]] += momentum_fluxes_negative_material[0] * one_cell_size;
               if constexpr( MF::IsEquationActive( Equation::MomentumY ) ) {
                  right_hand_side_positive_material[Equation::MomentumY][indices[0]][indices[1]][indices[2]] -= momentum_fluxes_positive_material[1] * one_cell_size;
                  right_hand_side_negative_material[Equation::MomentumY][indices[0]][indices[1]][indices[2]] += momentum_fluxes_negative_material[1] * one_cell_size;
               }
               if constexpr( MF::IsEquationActive( Equation::MomentumZ ) ) {
                  right_hand_side_positive_material[Equation::MomentumZ][indices[0]][indices[1]][indices[2]] -= momentum_fluxes_positive_material[2] * one_cell_size;
                  right_hand_side_negative_material[Equation::MomentumZ][indices[0]][indices[1]][indices[2]] += momentum_fluxes_negative_material[2] * one_cell_size;
               }
            }//if
         }   //k
      }      //j
   }         //i
}

/**
 * @brief      Adds the inviscid part to the interface stress tensor.
 *
 * @param      node                                    The considered node.
 * @param      interface_stress_tensor_positive_material  The interface stress tensor of the positive material.
 * @param      interface_stress_tensor_negative_material  The interface stress tensor of the negative material.
 */
void InterfaceStressTensorFluxes::AddInviscidPartToInterfaceStressTensor( Node const& node, double ( &interface_stress_tensor_positive_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )], double ( &interface_stress_tensor_negative_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const {

   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()]         = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
   double const( &interface_pressure_positive )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetInterfaceStateBuffer( InterfaceState::PressurePositive );
   double const( &interface_pressure_negative )[CC::TCX()][CC::TCY()][CC::TCZ()] = CC::CapillaryForcesActive() ?
                                                                                         node.GetInterfaceBlock().GetInterfaceStateBuffer( InterfaceState::PressureNegative ) :
                                                                                         node.GetInterfaceBlock().GetInterfaceStateBuffer( InterfaceState::PressurePositive );

   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            std::array<unsigned int, 3> const indices = { BIT::I2TX( i ), BIT::I2TY( j ), BIT::I2TZ( k ) };
            if( std::abs( interface_tags[indices[0]][indices[1]][indices[2]] ) <= ITTI( IT::NewCutCell ) ) {
               double const pressure_positive = interface_pressure_positive[indices[0]][indices[1]][indices[2]];
               double const pressure_negative = interface_pressure_negative[indices[0]][indices[1]][indices[2]];
               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  interface_stress_tensor_positive_material[i][j][k][r][r] -= pressure_positive;
                  interface_stress_tensor_negative_material[i][j][k][r][r] -= pressure_negative;
               }
            }//if
         }   //k
      }      //j
   }         //i
}

/**
 * @brief      Adds the viscous part to the interface stress tensor.
 * @param      node                                    The considered node.
 * @param      interface_stress_tensor_positive_material  The interface stress tensor of the positive material.
 * @param      interface_stress_tensor_negative_material  The interface stress tensor of the negative material.
 */
void InterfaceStressTensorFluxes::AddViscousPartToInterfaceStressTensor( Node const& node, double ( &interface_stress_tensor_positive_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )], double ( &interface_stress_tensor_negative_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const {

   double velocity_gradient[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )];
   double tau[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )];
   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
               for( unsigned int s = 0; s < DTI( CC::DIM() ); ++s ) {
                  velocity_gradient[i][j][k][r][s] = 0.0;
                  tau[i][j][k][r][s]               = 0.0;
               }// s
            }   // r
         }      // k
      }         // j
   }            // i

   CalculateVelocityGradientAtInterface( node, velocity_gradient );

   CalculateViscousStressTensor( node, velocity_gradient, tau );
   if constexpr( CC::Axisymmetric() ) {
      AddAxisymmetricPartToViscousStressTensor( node, tau );
   }

   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

   bool const is_negative_material_solid = material_manager_.IsSolidBoundary( negative_material_properties_.material_ );
   bool const is_positive_material_solid = material_manager_.IsSolidBoundary( positive_material_properties_.material_ );

   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            std::array<unsigned int, 3> const indices = { BIT::I2TX( i ), BIT::I2TY( j ), BIT::I2TZ( k ) };
            if( std::abs( interface_tags[indices[0]][indices[1]][indices[2]] ) <= ITTI( IT::NewCutCell ) ) {
               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  for( unsigned int s = 0; s < DTI( CC::DIM() ); ++s ) {
                     if( !is_negative_material_solid ) {
                        interface_stress_tensor_negative_material[i][j][k][r][s] += tau[i][j][k][r][s];
                     }
                     if( !is_positive_material_solid ) {
                        interface_stress_tensor_positive_material[i][j][k][r][s] += tau[i][j][k][r][s];
                     }
                  }// s
               }   // r
            }      //if
         }         //k
      }            //j
   }               //i
}

/**
 * @brief      Calculates the velocity gradient at the interface.
 * @param      node                            The node.
 * @param      velocity_gradient_at_interface  The velocity gradient at the interface as an indirect return parameter.
 */
void InterfaceStressTensorFluxes::CalculateVelocityGradientAtInterface( Node const& node, double ( &velocity_gradient_at_interface )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const {

   double const cell_size                                                = node.GetCellSize();
   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

   double real_material_velocity_x[CC::TCX()][CC::TCY()][CC::TCZ()];
   double real_material_velocity_y[CC::TCX()][CC::TCY()][CC::TCZ()];
   double real_material_velocity_z[CC::TCX()][CC::TCY()][CC::TCZ()];
   ComputeRealMaterialVelocity( node, real_material_velocity_x, real_material_velocity_y, real_material_velocity_z );

   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell ) ) {

               std::array<std::array<double, 3>, 3> const gradient = SU::JacobianMatrix<DerivativeStencilSetup::Concretize<viscous_fluxes_derivative_stencil_cell_center>::type>(
                     real_material_velocity_x, real_material_velocity_y, real_material_velocity_z, i, j, k, cell_size );

               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  for( unsigned int s = 0; s < DTI( CC::DIM() ); ++s ) {
                     velocity_gradient_at_interface[BIT::T2IX( i )][BIT::T2IY( j )][BIT::T2IZ( k )][r][s] = gradient[r][s];
                  }
               }
            }//if
         }   //k
      }      //j
   }         //i
}

/**
 * @brief      Calculates the viscous contribution to the stress tensor tau.
 * @param      node               The node.
 * @param[in]  velocity_gradient  The interface velocity gradient.
 * @param      tau                The viscous part of the stress tensor as an indirect return parameter.
 */
void InterfaceStressTensorFluxes::CalculateViscousStressTensor( Node const& node, double const ( &velocity_gradient )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )], double ( &tau )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const {

   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
   double const( &volume_fractions )[CC::TCX()][CC::TCY()][CC::TCZ()]    = node.GetInterfaceBlock().GetReinitializedBuffer( InterfaceDescription::VolumeFraction );

   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            std::array<unsigned int, 3> const indices = { BIT::I2TX( i ), BIT::I2TY( j ), BIT::I2TZ( k ) };
            if( std::abs( interface_tags[indices[0]][indices[1]][indices[2]] ) <= ITTI( IT::NewCutCell ) ) {
               std::array<double, 3> const interface_viscosity        = ComputeInterfaceViscosities( volume_fractions[indices[0]][indices[1]][indices[2]] );
               std::array<double, 3> const velocity_gradient_diagonal = { velocity_gradient[i][j][k][0][0], CC::DIM() != Dimension::One ? velocity_gradient[i][j][k][1][1] : 0.0, CC::DIM() == Dimension::Three ? velocity_gradient[i][j][k][2][2] : 0.0 };
               double const volume_viscosity_contribution             = interface_viscosity[2] * DimensionAwareConsistencyManagedSum( velocity_gradient_diagonal );

               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  for( unsigned int s = 0; s < DTI( CC::DIM() ); ++s ) {
                     tau[i][j][k][r][s] += interface_viscosity[0] * ( velocity_gradient[i][j][k][r][s] + velocity_gradient[i][j][k][s][r] );
                  }
               }

               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  tau[i][j][k][r][r] += volume_viscosity_contribution;
               }
            }//if
         }   //k
      }      //j
   }         //i
}

/**
 * @brief      Adds the axisymmetric viscous part to the stress tensor
 *
 * @param      node               The node.
 * @param      tau                The viscid part of the stress tensor as an indirect return parameter.
 */
void InterfaceStressTensorFluxes::AddAxisymmetricPartToViscousStressTensor( Node const& node, double ( &tau )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const {

   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
   double const( &volume_fractions )[CC::TCX()][CC::TCY()][CC::TCZ()]    = node.GetInterfaceBlock().GetReinitializedBuffer( InterfaceDescription::VolumeFraction );

   double const cell_size = node.GetCellSize();

   double real_fluid_velocity_x[CC::TCX()][CC::TCY()][CC::TCZ()];
   double real_fluid_velocity_y[CC::TCX()][CC::TCY()][CC::TCZ()];
   double real_fluid_velocity_z[CC::TCX()][CC::TCY()][CC::TCZ()];
   ComputeRealMaterialVelocity( node, real_fluid_velocity_x, real_fluid_velocity_y, real_fluid_velocity_z );

   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            std::array<unsigned int, 3> const indices = { BIT::I2TX( i ), BIT::I2TY( j ), BIT::I2TZ( k ) };
            if( std::abs( interface_tags[indices[0]][indices[1]][indices[2]] ) <= ITTI( IT::NewCutCell ) ) {
               std::array<double, 3> const interface_viscosity = ComputeInterfaceViscosities( volume_fractions[indices[0]][indices[1]][indices[2]] );

               double const radius                        = std::get<0>( node.GetBlockCoordinates() ) + ( static_cast<double>( i ) + 0.5 ) * cell_size;
               double const volume_viscosity_contribution = interface_viscosity[2] * real_fluid_velocity_x[indices[0]][indices[1]][indices[2]] / radius;

               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  tau[i][j][k][r][r] += volume_viscosity_contribution;
               }
            }//if
         }   //k
      }      //j
   }         //i
}

/**
 * @brief      Calculates the real-material velocity field.
 *
 * @param      node                   The node.
 * @param      real_material_velocity_x  The real-material velocity field in x direction.
 * @param      real_material_velocity_y  The real-material velocity field in y direction.
 * @param      real_material_velocity_z  The real-material velocity field in z direction.
 */
void InterfaceStressTensorFluxes::ComputeRealMaterialVelocity( Node const& node, double ( &real_material_velocity_x )[CC::TCX()][CC::TCY()][CC::TCZ()], double ( &real_material_velocity_y )[CC::TCX()][CC::TCY()][CC::TCZ()], double ( &real_material_velocity_z )[CC::TCX()][CC::TCY()][CC::TCZ()] ) const {

   double const( &volume_fraction )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer( InterfaceDescription::VolumeFraction );

   Block const& positive_material = node.GetPhaseByMaterial( positive_material_properties_.material_ );
   Block const& negative_material = node.GetPhaseByMaterial( negative_material_properties_.material_ );

   PrimeStates const& positive_states = positive_material.GetPrimeStateBuffer();
   PrimeStates const& negative_states = negative_material.GetPrimeStateBuffer();

   for( unsigned int i = 0; i < CC::TCX(); ++i ) {
      for( unsigned int j = 0; j < CC::TCY(); ++j ) {
         for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
            real_material_velocity_x[i][j][k] = volume_fraction[i][j][k] * positive_states[PrimeState::VelocityX][i][j][k] + ( 1.0 - volume_fraction[i][j][k] ) * negative_states[PrimeState::VelocityX][i][j][k];
            real_material_velocity_y[i][j][k] = CC::DIM() != Dimension::One ? volume_fraction[i][j][k] * positive_states[PrimeState::VelocityY][i][j][k] + ( 1.0 - volume_fraction[i][j][k] ) * negative_states[PrimeState::VelocityY][i][j][k] : 0.0;
            real_material_velocity_z[i][j][k] = CC::DIM() == Dimension::Three ? volume_fraction[i][j][k] * positive_states[PrimeState::VelocityZ][i][j][k] + ( 1.0 - volume_fraction[i][j][k] ) * negative_states[PrimeState::VelocityZ][i][j][k] : 0.0;
         }
      }
   }
}
