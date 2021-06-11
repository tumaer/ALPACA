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
#include <catch2/catch.hpp>

#include <vector>
#include <algorithm>

#include "materials/equations_of_state/generic_stiffened_gas.h"

SCENARIO( "Pressure Computation", "[1rank]" ) {

   GIVEN( "An unsafe set of variables" ) {
      constexpr double density    = 0.0;
      constexpr double momentum_x = 0.1;
      constexpr double momentum_y = 0.2;
      constexpr double momentum_z = 0.3;
      constexpr double energy     = 2.0;
      constexpr double gamma      = 1.4;
      constexpr double pi         = 0.0;
      WHEN( "We compute the speed of sound in a safe way" ) {
         constexpr bool safely = true;
         THEN( "The result is positive." ) {
            REQUIRE( GenericStiffenedGas::CalculatePressure<safely>( density, momentum_x, momentum_y, momentum_z, energy, gamma, pi ) > 0.0 );
         }
      }
      WHEN( "We compute the pressure in an unsafe way" ) {
         constexpr bool safely = false;
         THEN( "No proper number is obtained." ) {
            REQUIRE_FALSE( std::isnormal( GenericStiffenedGas::CalculatePressure<safely>( density, momentum_x, momentum_y, momentum_z, energy, gamma, pi ) ) );
         }
      }
   }
}

SCENARIO( "Energy Computation", "[1rank]" ) {

   GIVEN( "An arbitrary set of variables" ) {
      constexpr double density    = 1.1;
      constexpr double velocity_x = 0.1;
      constexpr double velocity_y = 0.2;
      constexpr double velocity_z = 0.3;
      constexpr double pressure   = 2.0;
      constexpr double gamma      = 1.4;
      constexpr double pi         = 3.0;
      WHEN( "We compute the energy for a stiffened equation of state" ) {
         THEN( "The computed energy matches the expectation" ) {
            REQUIRE( GenericStiffenedGas::CalculateEnergy( density, velocity_x, velocity_y, velocity_z, pressure, gamma, pi ) == Approx( 15.577 ).margin( 1.0e-15 ) );
         }
      }
   }
}

SCENARIO( "Speed of Sound Computation", "[1rank]" ) {

   GIVEN( "An unsafe set of variables" ) {
      constexpr double density  = 1.0;
      constexpr double pressure = -0.0001;
      constexpr double gamma    = 1.4;
      constexpr double pi       = 0.0;
      WHEN( "We compute the speed of sound in a safe way" ) {
         constexpr bool safely = true;
         THEN( "The result is positive." ) {
            REQUIRE( GenericStiffenedGas::CalculateSpeedOfSound<safely>( density, pressure, gamma, pi ) > 0.0 );
         }
      }
      WHEN( "We compute the speed of sound in an unsafe way" ) {
         constexpr bool safely = false;
         THEN( "A NaN is triggered." ) {
            REQUIRE( std::isnan( GenericStiffenedGas::CalculateSpeedOfSound<safely>( density, pressure, gamma, pi ) ) );
         }
      }
   }
}
