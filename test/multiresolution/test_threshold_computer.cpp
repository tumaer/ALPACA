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
#include "multiresolution/threshold_computer.h"
#include "enums/dimension_definition.h"

SCENARIO( "Multiresolution thresholds are correctly computed", "[1rank]" ) {

   GIVEN( "1D setup and a convergence order of zero" ) {
      constexpr Dimension dim                  = Dimension::One;
      constexpr unsigned int convergence_order = 0;
      WHEN( "The maximum level is equal the reference level" ) {
         constexpr unsigned int maximum_level   = 6;
         constexpr unsigned int reference_level = maximum_level;
         THEN( "The thresholds follows the equation eps_ref * 2 ^ ( lmax - l ) " ) {
            constexpr double user_epsilon_ref                     = 1.0;
            constexpr double thresholder_epsilon_ref              = 1.0;// Thresholder computed epsilon_ref (derived from user value)
            ThresholdComputer<dim, convergence_order> thresholder = ThresholdComputer<dim, convergence_order>( maximum_level, reference_level, user_epsilon_ref );
            REQUIRE( thresholder.ThresholdOnLevel( 0 ) == Approx( std::pow( 2.0, -6.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 1 ) == Approx( std::pow( 2.0, -5.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 2 ) == Approx( std::pow( 2.0, -4.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 3 ) == Approx( std::pow( 2.0, -3.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 4 ) == Approx( std::pow( 2.0, -2.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 5 ) == Approx( std::pow( 2.0, -1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 6 ) == Approx( std::pow( 2.0, 0.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 7 ) == Approx( std::pow( 2.0, 1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 8 ) == Approx( std::pow( 2.0, 2.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 9 ) == Approx( std::pow( 2.0, 3.0 ) * thresholder_epsilon_ref ) );
         }
      }
      WHEN( "The maximum level is larger than the reference level" ) {
         constexpr unsigned int maximum_level   = 6;
         constexpr unsigned int reference_level = 4;
         THEN( "The thresholds follows the equation eps_ref * 2 ^ - ( lmax - lref )  * 2 ^ - ( lmax - l ) " ) {
            constexpr double user_epsilon_ref                     = 0.5;
            constexpr double thresholder_epsilon_ref              = 0.125;// Thresholder computed epsilon_ref (derived from user value)
            ThresholdComputer<dim, convergence_order> thresholder = ThresholdComputer<dim, convergence_order>( maximum_level, reference_level, user_epsilon_ref );
            // 0.125 should be the internal used epsilon_ref != provided epsilon_ref
            REQUIRE( thresholder.ThresholdOnLevel( 0 ) == Approx( std::pow( 2.0, -6.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 1 ) == Approx( std::pow( 2.0, -5.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 2 ) == Approx( std::pow( 2.0, -4.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 3 ) == Approx( std::pow( 2.0, -3.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 4 ) == Approx( std::pow( 2.0, -2.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 5 ) == Approx( std::pow( 2.0, -1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 6 ) == Approx( std::pow( 2.0, 0.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 7 ) == Approx( std::pow( 2.0, 1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 8 ) == Approx( std::pow( 2.0, 2.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 9 ) == Approx( std::pow( 2.0, 3.0 ) * thresholder_epsilon_ref ) );
         }
      }
      WHEN( "The maximum level is smaller than the reference level" ) {
         constexpr unsigned int maximum_level   = 6;
         constexpr unsigned int reference_level = 8;
         THEN( "The thresholds follows the equation eps_ref * 2 ^ - ( lmax - lref )  * 2 ^ - ( lmax - l ) " ) {
            constexpr double user_epsilon_ref                     = 1.0;
            constexpr double thresholder_epsilon_ref              = 4.0;// Thresholder computed epsilon_ref (derived from user value)
            ThresholdComputer<dim, convergence_order> thresholder = ThresholdComputer<dim, convergence_order>( maximum_level, reference_level, user_epsilon_ref );
            // 4.0 should be the internal used epsilon_ref != provided epsilon_ref
            REQUIRE( thresholder.ThresholdOnLevel( 0 ) == Approx( std::pow( 2.0, -6.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 1 ) == Approx( std::pow( 2.0, -5.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 2 ) == Approx( std::pow( 2.0, -4.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 3 ) == Approx( std::pow( 2.0, -3.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 4 ) == Approx( std::pow( 2.0, -2.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 5 ) == Approx( std::pow( 2.0, -1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 6 ) == Approx( std::pow( 2.0, 0.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 7 ) == Approx( std::pow( 2.0, 1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 8 ) == Approx( std::pow( 2.0, 2.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 9 ) == Approx( std::pow( 2.0, 3.0 ) * thresholder_epsilon_ref ) );
         }
      }
   }

   GIVEN( "The maximum level is two larger than the reference level" ) {
      constexpr unsigned int reference_level = 5;
      constexpr unsigned int maximum_level   = reference_level + 2;
      WHEN( "The convergence order is two and the dimension is two" ) {
         constexpr unsigned int convergence_order = 2;
         constexpr Dimension dim                  = Dimension::Two;
         THEN( "The thresholds follows the equation eps_ref * 2 ^ - 3 * ( lmax - lref )  * 2 ^ - 2 * ( lmax - l ) " ) {
            constexpr double user_epsilon_ref                     = 1.0;
            constexpr double thresholder_epsilon_ref              = 0.015625;// Thresholder computed epsilon_ref (derived from user value)
            ThresholdComputer<dim, convergence_order> thresholder = ThresholdComputer<dim, convergence_order>( maximum_level, reference_level, user_epsilon_ref );
            // 0.015625 should be the internal used epsilon_ref != provided epsilon_ref
            REQUIRE( thresholder.ThresholdOnLevel( 0 ) == Approx( std::pow( 2.0, 2.0 * -7.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 1 ) == Approx( std::pow( 2.0, 2.0 * -6.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 2 ) == Approx( std::pow( 2.0, 2.0 * -5.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 3 ) == Approx( std::pow( 2.0, 2.0 * -4.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 4 ) == Approx( std::pow( 2.0, 2.0 * -3.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 5 ) == Approx( std::pow( 2.0, 2.0 * -2.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 6 ) == Approx( std::pow( 2.0, 2.0 * -1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 7 ) == Approx( std::pow( 2.0, 2.0 * 0.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 8 ) == Approx( std::pow( 2.0, 2.0 * 1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 9 ) == Approx( std::pow( 2.0, 2.0 * 2.0 ) * thresholder_epsilon_ref ) );
         }
      }

      WHEN( "The convergence order is five and the dimension is three" ) {
         constexpr unsigned int convergence_order = 5;
         constexpr Dimension dim                  = Dimension::Three;
         THEN( "The thresholds follows the equation eps_ref * 2 ^ - 6 * ( lmax - lref )  * 2 ^ - 3 * ( lmax - l ) " ) {
            constexpr double user_epsilon_ref                     = 1.0;
            constexpr double thresholder_epsilon_ref              = 0.000244140625;// Thresholder computed epsilon_ref (derived from user value)
            ThresholdComputer<dim, convergence_order> thresholder = ThresholdComputer<dim, convergence_order>( maximum_level, reference_level, user_epsilon_ref );

            REQUIRE( thresholder.ThresholdOnLevel( 0 ) == Approx( std::pow( 2.0, 3.0 * -7.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 1 ) == Approx( std::pow( 2.0, 3.0 * -6.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 2 ) == Approx( std::pow( 2.0, 3.0 * -5.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 3 ) == Approx( std::pow( 2.0, 3.0 * -4.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 4 ) == Approx( std::pow( 2.0, 3.0 * -3.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 5 ) == Approx( std::pow( 2.0, 3.0 * -2.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 6 ) == Approx( std::pow( 2.0, 3.0 * -1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 7 ) == Approx( std::pow( 2.0, 3.0 * 0.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 8 ) == Approx( std::pow( 2.0, 3.0 * 1.0 ) * thresholder_epsilon_ref ) );
            REQUIRE( thresholder.ThresholdOnLevel( 9 ) == Approx( std::pow( 2.0, 3.0 * 2.0 ) * thresholder_epsilon_ref ) );
         }
      }
   }
}
