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
#ifndef INTERFACE_STRESS_TENSOR_FLUXES_H
#define INTERFACE_STRESS_TENSOR_FLUXES_H

#include "materials/equation_of_state_definitions.h"
#include "materials/material_manager.h"
#include "topology/node.h"
#include <vector>

/**
 * @brief Class for interface stress tensor fluxes. Calculates the interface fluxes caused by the stress tensor. Thus, it considers inviscid and viscous contributions.
 */
class InterfaceStressTensorFluxes {

   /**
    * @brief A struct to bundle viscous flux related information of a material.
    */
   struct ViscousMaterialProperties {
      MaterialName const material_;
      double const mu_shear_;
      double const mu_bulk_;

      /**
       * @brief      Default constructor of the struct.
       *
       * @param[in]  material  The material for which viscosities are saved.
       * @param[in]  mu        A vector containing the shear and bulk viscosity of the material.
       */
      ViscousMaterialProperties( MaterialName const material, std::vector<double> const mu ) : material_( material ),
                                                                                               mu_shear_( mu[0] ),
                                                                                               mu_bulk_( mu[1] ) { /* Empty, besides initializer list */
      }
   };

   MaterialManager const& material_manager_;
   ViscousMaterialProperties const positive_material_properties_;
   ViscousMaterialProperties const negative_material_properties_;

   static constexpr double epsilon_ = std::numeric_limits<double>::epsilon();

   std::array<double, 3> ComputeInterfaceViscosities( double const volume_fraction ) const;

   void AddInviscidPartToInterfaceStressTensor( Node const& node, double ( &interface_stress_tensor_positive_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )], double ( &interface_stress_tensor_negative_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const;
   void AddViscousPartToInterfaceStressTensor( Node const& node, double ( &interface_stress_tensor_positive_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )], double ( &interface_stress_tensor_negative_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const;

   void CalculateVelocityGradientAtInterface( Node const& node, double ( &velocity_gradient_at_interface )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const;

   void CalculateViscousStressTensor( Node const& node, double const ( &velocity_gradient )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )], double ( &tau )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const;

   void AddAxisymmetricPartToViscousStressTensor( Node const& node, double ( &tau )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const;

   void AddFluxesToRightHandSide( Node& node, double const ( &delta_aperture_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3], double const ( &u_interface_normal_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3], double const ( &interface_stress_tensor_positive_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )], double const ( &interface_stress_tensor_negative_material )[CC::ICX()][CC::ICY()][CC::ICZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) const;

   void ComputeRealMaterialVelocity( Node const& node, double ( &real_material_velocity_x )[CC::TCX()][CC::TCY()][CC::TCZ()], double ( &real_material_velocity_y )[CC::TCX()][CC::TCY()][CC::TCZ()], double ( &real_material_velocity_z )[CC::TCX()][CC::TCY()][CC::TCZ()] ) const;

public:
   InterfaceStressTensorFluxes() = delete;
   explicit InterfaceStressTensorFluxes( MaterialManager const& material_manager,
                                         MaterialName const material_positive, std::vector<double> const mu_positive,
                                         MaterialName const material_negative, std::vector<double> const mu_negative );
   ~InterfaceStressTensorFluxes()                                    = default;
   InterfaceStressTensorFluxes( InterfaceStressTensorFluxes const& ) = delete;
   InterfaceStressTensorFluxes& operator=( InterfaceStressTensorFluxes const& ) = delete;
   InterfaceStressTensorFluxes( InterfaceStressTensorFluxes&& )                 = delete;
   InterfaceStressTensorFluxes& operator=( InterfaceStressTensorFluxes&& ) = delete;

   void ComputeInterfaceFluxes( Node& node, double const ( &delta_aperture_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3], double const ( &u_interface_normal_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3] ) const;
};

#endif//INTERFACE_STRESS_TENSOR_FLUXES_H
