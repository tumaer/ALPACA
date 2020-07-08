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
#include <catch.hpp>
#include "solvers/source_term_contributions/axisymmetric_viscous_volume_forces.h"
#include "topology/node.h"
#include "materials/equations_of_state/stiffened_gas.h"

namespace {
   MaterialManager ReturnMaterialManagerWithViscosities( double const shear_viscosity, double const bulk_viscosity ) {
      // Initialize the unit handler class
      UnitHandler const unit_handler( 1.0, 1.0, 1.0, 1.0 );

      // Here the stiffened gas complete safe equation of state is used since it provides a temperature computation and the unit test does not depend on the
      // activation of the temperature in the primestate struct
      std::unordered_map<std::string, double> const eos_data = { { "gamma", 1.4 }, { "backgroundPressure", 0.0 } };
      std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGas const>( eos_data, unit_handler ) );

      // Define material properties and initialize material
      double const specific_heat_capacity = 3.0;
      double const thermal_heat_conductivity = 4.0;

      // Instantiate material
      std::vector<Material> materials;
      materials.emplace_back( Material( std::move( equation_of_state ), bulk_viscosity, shear_viscosity, thermal_heat_conductivity, specific_heat_capacity,
                                       nullptr, nullptr, unit_handler ) );

      // Instantiate material pairing
      std::vector<MaterialPairing> material_pairings;

      return MaterialManager( std::move( materials ), std::move( material_pairings ) );
   }
}

SCENARIO( "Volume Forces Calculation Correctness", "[1rank]" ) {

   GIVEN( "An axisymmetric viscous volume forces calculator" ) {
      std::pair<MaterialName const, Block> mat_block( std::piecewise_construct, std::make_tuple( MaterialName::MaterialOne ), std::make_tuple() );
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               mat_block.second.GetPrimeStateBuffer( PrimeState::VelocityX )[i][j][k] = static_cast<double>( i ) * 3.0;
               mat_block.second.GetPrimeStateBuffer( PrimeState::VelocityY )[i][j][k] = static_cast<double>( j ) * 4.0;
            }
         }
      }

      constexpr double cell_size = 1.0;
      WHEN( "The shear and bulk viscosity are zero" ) {
         constexpr double node_origin_x = 0.0;
         MaterialManager const material_manager = ReturnMaterialManagerWithViscosities( 0.0, 0.0 );
         AxisymmetricViscousVolumeForces forces_calculator = AxisymmetricViscousVolumeForces( material_manager );
         double volume_forces[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                     volume_forces[e][i][j][k] = 0.0;
                  }
               } //j
            } //i
         } //equation

         forces_calculator.ComputeForces( mat_block, volume_forces, cell_size, node_origin_x );

         THEN( "The volume forces are zero" ) {
             for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
                for( unsigned int i = 0; i < CC::ICX(); ++i ) {
                   for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                      REQUIRE( volume_forces[e][i][j][0] == Approx( 0.0 ) );
                   } //j
                } //i
             } //equation
         }
      }
      WHEN( "Viscosities are twice different" ) {
         double const node_origin_x = 0.0;
         double first_volume_forces[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         double second_volume_forces[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                     first_volume_forces[e][i][j][k] = 0.0;
                     second_volume_forces[e][i][j][k] = 0.0;
                  }
               } //j
            } //i
         } //equation
         MaterialManager const first_material_manager = ReturnMaterialManagerWithViscosities( 3.4, 0.0 );
         AxisymmetricViscousVolumeForces first_forces_calculator = AxisymmetricViscousVolumeForces( first_material_manager );
         first_forces_calculator.ComputeForces(mat_block, first_volume_forces, cell_size, node_origin_x);
         MaterialManager const second_material_manager = ReturnMaterialManagerWithViscosities( 6.8, 0.0 );
         AxisymmetricViscousVolumeForces second_forces_calculator = AxisymmetricViscousVolumeForces( second_material_manager );
         second_forces_calculator.ComputeForces(mat_block, second_volume_forces, cell_size, node_origin_x);
         THEN( "Volume forces are twice different" ) {
             for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
                for( unsigned int i = 0; i < CC::ICX(); ++i ) {
                   for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                       REQUIRE( 2.0 * first_volume_forces[e][i][j][0] == Approx( second_volume_forces[e][i][j][0] ) );
                   } //j
                } //i
             } //equation
         }
      }
      WHEN( "Cell-center coordinates are twice different due to the difference of the block position" ) {
         constexpr double first_node_origin_x = 0.0;
         constexpr double second_node_origin_x = 4.5;
         double first_volume_forces[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         double second_volume_forces[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];
         for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                     first_volume_forces[e][i][j][k] = 0.0;
                     second_volume_forces[e][i][j][k] = 0.0;
                  }
               } //j
            } //i
         } //equation
         MaterialManager const material_manager = ReturnMaterialManagerWithViscosities( 5.0, 0.0 );
         AxisymmetricViscousVolumeForces forces_calculator = AxisymmetricViscousVolumeForces( material_manager );
         forces_calculator.ComputeForces( mat_block, first_volume_forces, cell_size, first_node_origin_x );
         forces_calculator.ComputeForces( mat_block, second_volume_forces, cell_size, second_node_origin_x );

         THEN( "r-momentum forces are twice different" ) {
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  REQUIRE( first_volume_forces[ETI( Equation::MomentumY )][i][j][0] == Approx( 2.0 * second_volume_forces[ETI( Equation::MomentumY )][i][j][0] ) );
               } //j
            } //i
         }
      }
   }
}