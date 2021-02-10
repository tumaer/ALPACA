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

#include <vector>
#include <string>
#include <tuple>

#include "initial_condition/prime_state_initializer.h"
#include "unit_handler.h"
#include "initial_condition/test_initial_condition_helper.h"

namespace {

   std::tuple<std::vector<std::string>, std::vector<std::string>> CreateMaterialInitialConditions( std::vector<MaterialName> const materials,
                                                                                                   std::array<double, 3> const factors ) {
      std::vector<std::string> material_initial_expressions( materials.size(), "" );
      std::vector<std::string> prime_state_variable_names( MF::ANOP() );
      for( PrimeState const p : MF::ASOP() ) {
         std::string const prime_state_name   = StringOperations::RemoveSpaces( std::string( MF::InputName( p ) ) );
         prime_state_variable_names[PTI( p )] = prime_state_name;
         if( !prime_state_name.empty() ) {
            for( auto const material : materials ) {
               material_initial_expressions[MTI( material )] += prime_state_name + " := " + std::to_string( double( MTI( material ) + 1 ) * double( PTI( p ) + 1 ) ) + "(" +
                                                                std::to_string( factors[0] ) + " * x + " +
                                                                std::to_string( factors[1] ) + " * y + " +
                                                                std::to_string( factors[2] ) + " * z);";
            }
         }
      }

      return std::make_tuple( prime_state_variable_names, material_initial_expressions );
   }
}// namespace

SCENARIO( "The prime state initializer works properly.", "[1rank]" ) {
   GIVEN( "A topology with: node_size = 1.0, maximum level = 0, n_materials = 2. The prime state expression are linear in each direction and differ"
          " in a single factor for each prime state and material. And a unit handler." ) {
      // General data
      constexpr double node_size                        = 1.0;
      constexpr unsigned int max_level                  = 0;
      std::vector<std::vector<nid_t>> const level_nodes = TestInitialCondition::GetChildNodesForSingleRootNode( max_level );
      std::vector<MaterialName> const materials         = { MaterialName::MaterialOne, MaterialName::MaterialTwo };

      WHEN( "The unit handler has factors of unity." ) {
         // Create the unit handler
         UnitHandler const unit_handler( 1.0, 1.0, 1.0, 1.0 );
         // Create the initial conditions
         std::tuple<std::vector<std::string>, std::vector<std::string>> const primes_and_initial_conditions( CreateMaterialInitialConditions( materials, { 0.5, 1.0, 1.5 } ) );
         // Create the initializer
         PrimeStateInitializer initializer( PrimeStateInitializer( std::get<1>( primes_and_initial_conditions ),
                                                                   std::get<0>( primes_and_initial_conditions ),
                                                                   node_size,
                                                                   unit_handler ) );
         // Compute the primes for each material
         double buffer_mat1[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()];
         double buffer_mat2[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()];
         initializer.GetInitialPrimeStates( level_nodes[0][0], MaterialName::MaterialOne, buffer_mat1 );
         initializer.GetInitialPrimeStates( level_nodes[0][0], MaterialName::MaterialTwo, buffer_mat2 );

         THEN( "The prime state values of the two materials must differ in the material index as the factor." ) {
            double const conversion_factor = double( MTI( MaterialName::MaterialTwo ) + 1 ) / double( MTI( MaterialName::MaterialOne ) + 1 );
            // Loop through all elements and check the data
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                     for( PrimeState const p : MF::ASOP() ) {
                        REQUIRE( ( conversion_factor * buffer_mat1[PTI( p )][i][j][k] ) == buffer_mat2[PTI( p )][i][j][k] );
                     }
                  }
               }
            }
         }
         THEN( "The prime states for a single material are the same divided by its prime state index, but not zero. Zero if no name exists for the prime state." ) {
            // Loop through all elements and check the data
            for( unsigned int i = 0; i < CC::ICX(); ++i ) {
               for( unsigned int j = 0; j < CC::ICY(); ++j ) {
                  for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                     for( PrimeState const p : MF::ASOP() ) {
                        if( std::string( MF::InputName( p ) ).empty() ) {
                           REQUIRE( buffer_mat1[PTI( p )][i][j][k] == 0.0 );
                           REQUIRE( buffer_mat2[PTI( p )][i][j][k] == 0.0 );
                        } else {
                           double const factor = double( PTI( PrimeState::Density ) + 1 ) / double( PTI( p ) + 1 );
                           REQUIRE( buffer_mat1[PTI( PrimeState::Density )][i][j][k] == factor * buffer_mat1[PTI( p )][i][j][k] );
                           REQUIRE( buffer_mat1[PTI( p )][i][j][k] != 0.0 );
                           REQUIRE( buffer_mat2[PTI( PrimeState::Density )][i][j][k] == factor * buffer_mat2[PTI( p )][i][j][k] );
                           REQUIRE( buffer_mat2[PTI( p )][i][j][k] != 0.0 );
                        }
                     }
                  }
               }
            }
         }
      }
   }
}
