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

#include "user_expression.h"
#include "utilities/vector_utilities.h"

SCENARIO( "User expression evaluation", "[1rank]" ) {

   GIVEN( "An arbitrary user expression." ) {
      WHEN( "The point contains less elements then the function input arguments." ) {
         std::vector<double> point( 1, 0.0 );
         std::vector<std::string> const variables_out = { "f" };
         std::vector<std::string> const variables_in  = { "x", "y" };
         std::string const expression                 = "f := x* y";
         THEN( "An exception should be thrown in the constructor of the UserExpression." ) {
            REQUIRE_THROWS( UserExpression( expression, variables_out, variables_in, point ) );
         }
      }
      WHEN( "The function input arguments contains less elements then point contains." ) {
         std::vector<double> point( 2, 0.0 );
         std::vector<std::string> const variables_out = { "f" };
         std::vector<std::string> const variables_in  = { "x" };
         std::string const expression                 = "f := x";
         THEN( "An exception should be thrown in the constructor of the UserExpression." ) {
            REQUIRE_THROWS( UserExpression( expression, variables_out, variables_in, point ) );
         }
      }
   }

   GIVEN( "A linear user expression of type f=2*x+5 for single corodinate computation." ) {
      std::vector<double> point( 1, 0.0 );
      std::vector<std::string> const variables_out = { "f" };
      std::vector<std::string> const variables_in  = { "x" };
      std::string const expression                 = "f := 2 * x + 5";
      UserExpression user_expression               = UserExpression( expression, variables_out, variables_in, point );

      WHEN( "The point is {x} = {0.0}" ) {
         point[0] = 0.0;
         THEN( "f should be approximately 5.0" ) {
            REQUIRE( user_expression.GetValue( "f" ) == Approx( 5.0 ) );
         }
      }
      WHEN( "The point is {x} = {-2.5}" ) {
         point[0] = -2.5;
         THEN( "f should be approximately 0.0" ) {
            REQUIRE( user_expression.GetValue( "f" ) == Approx( 0.0 ) );
         }
      }
      WHEN( "The point is {x} = {-2.625}" ) {
         point[0] = -2.625;
         THEN( "f should be approximately -0.25" ) {
            REQUIRE( user_expression.GetValue( "f" ) == Approx( -0.25 ) );
         }
      }
   }
   GIVEN( "A quadratic user expression of type f=x*x+2.5, g=x*y for multiple corodinate computation." ) {
      std::vector<double> point( 2, 0.0 );
      std::vector<std::string> const variables_out = { "f", "g" };
      std::vector<std::string> const variables_in  = { "x", "y" };
      std::string const expression                 = "f := x * x + 2.5; g := x * y;";
      UserExpression user_expression               = UserExpression( expression, variables_out, variables_in, point );

      WHEN( "The point is {x,y} = {0.0, 0.0}" ) {
         point[0] = 0.0;
         point[1] = 0.0;
         THEN( "Then, approximately should hold {f,g} = {2.5, 0.0}" ) {
            REQUIRE( user_expression.GetValue( "f" ) == Approx( 2.5 ) );
            REQUIRE( user_expression.GetValue( "g" ) == Approx( 0.0 ) );
         }
      }
      WHEN( "The point is {x,y} = {0.0, 2.125}" ) {
         point[0] = 0.0;
         point[1] = 2.125;
         THEN( "Then, approximately should hold {f,g} = {2.5, 0.0}" ) {
            REQUIRE( user_expression.GetValue( "f" ) == Approx( 2.5 ) );
            REQUIRE( user_expression.GetValue( "g" ) == Approx( 0.0 ) );
         }
      }
      WHEN( "The point is {x,y} = {1.5, 0.0}" ) {
         point[0] = 1.5;
         point[1] = 0.0;
         THEN( "Then, approximately should hold {f,g} = {4.75, 0.0}" ) {
            REQUIRE( user_expression.GetValue( "f" ) == Approx( 4.75 ) );
            REQUIRE( user_expression.GetValue( "g" ) == Approx( 0.0 ) );
         }
      }
      WHEN( "The point is {x,y} = {2.75, -1.25}" ) {
         point[0] = 2.75;
         point[1] = -1.25;
         THEN( "Then, approximately should hold {f,g} = {10.0625, -3.4375}" ) {
            REQUIRE( user_expression.GetValue( "f" ) == Approx( 10.0625 ) );
            REQUIRE( user_expression.GetValue( "g" ) == Approx( -3.4375 ) );
         }
      }
   }
   GIVEN( "A parametric user expression of type x=r*sin(theta)*cos(phi), y=r*sin(theta)*sin(phi), z=r*cos(theta) for spherical corodinate computation." ) {
      std::vector<double> point( 3, 0.0 );
      std::vector<std::string> const variables_out = { "x", "y", "z" };
      std::vector<std::string> const variables_in  = { "r", "theta", "phi" };
      std::string const expression                 = "x := r*sin(theta)*cos(phi); y := r*sin(theta)*sin(phi); z := r*cos(theta);";
      UserExpression user_expression               = UserExpression( expression, variables_out, variables_in, point );

      WHEN( "The radius r has three arbitrary positive values and the angles have 5 points in the range theta = [-pi, pi) and phi = [0, 2pi)" ) {
         constexpr std::array<double, 3> radius = { 0.5, 2.25, 5.5 };
         constexpr double delta_angle           = 2 * 3.1415926 / 5;
         std::array<std::vector<std::array<double, 3>>, 3> coordinates;
         for( unsigned int r = 0; r < 3; r++ ) {
            point[0] = radius[r];
            for( unsigned int n_theta = 0; n_theta < 5; n_theta++ ) {
               point[1] = 0.0 + n_theta * delta_angle;
               for( unsigned int n_phi = 0; n_phi < 5; n_phi++ ) {
                  point[2] = -3.1415926 + n_phi * delta_angle;
                  coordinates[r].push_back( { user_expression.GetValue( "x" ), user_expression.GetValue( "y" ), user_expression.GetValue( "z" ) } );
               }
            }
         }
         THEN( "The coordinates distance between two radii should be same for all angles." ) {
            double const dist_r1_r2 = VU::Distance( coordinates[0][0], coordinates[1][0] );
            for( unsigned int i = 1; i < coordinates[0].size(); i++ ) {
               REQUIRE( VU::Distance( coordinates[0][i], coordinates[1][i] ) == Approx( dist_r1_r2 ) );
            }
            double const dist_r1_r3 = VU::Distance( coordinates[0][0], coordinates[2][0] );
            for( unsigned int i = 1; i < coordinates[0].size(); i++ ) {
               REQUIRE( VU::Distance( coordinates[0][i], coordinates[2][i] ) == Approx( dist_r1_r3 ) );
            }
            double const dist_r2_r3 = VU::Distance( coordinates[1][0], coordinates[2][0] );
            for( unsigned int i = 1; i < coordinates[1].size(); i++ ) {
               REQUIRE( VU::Distance( coordinates[1][i], coordinates[2][i] ) == Approx( dist_r2_r3 ) );
            }
         }
         THEN( "The coordinates distance for all points between the point and center at {0,0,0} must be the same for each radius." ) {
            double const dist_r1 = VU::Distance( coordinates[0][0], { 0.0, 0.0, 0.0 } );
            for( unsigned int i = 1; i < coordinates[0].size(); i++ ) {
               REQUIRE( VU::Distance( coordinates[0][i], { 0.0, 0.0, 0.0 } ) == Approx( dist_r1 ) );
            }
            double const dist_r2 = VU::Distance( coordinates[1][0], { 0.0, 0.0, 0.0 } );
            for( unsigned int i = 1; i < coordinates[0].size(); i++ ) {
               REQUIRE( VU::Distance( coordinates[1][i], { 0.0, 0.0, 0.0 } ) == Approx( dist_r2 ) );
            }
            double const dist_r3 = VU::Distance( coordinates[2][0], { 0.0, 0.0, 0.0 } );
            for( unsigned int i = 1; i < coordinates[0].size(); i++ ) {
               REQUIRE( VU::Distance( coordinates[2][i], { 0.0, 0.0, 0.0 } ) == Approx( dist_r3 ) );
            }
         }
      }
   }
}
