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

#include <array>
#include <tinyxml2.h>

#include "input_output/utilities/stl_utilities.h"
#include "utilities/vector_utilities.h"

SCENARIO( "The triangle creation works properly", "[1rank]" ) {
   GIVEN( "An equilateral triangle with fixed non-zero z-coordinate and edge length one." ) {
      constexpr std::array<double, 3> p1 = { 1.0, 1.5, 2.0 };
      constexpr std::array<double, 3> p2 = { 1.5, 2.366025404, 2.0 };
      constexpr std::array<double, 3> p3 = { 2.0, 1.5, 2.0 };
      WHEN( "The normal vector points in positive z-direction and is not normalized." ) {
         constexpr std::array<double, 3> normal = { 0.0, 0.0, 2.34567 };
         StlUtilities::Triangle const triangle( normal, p1, p2, p3 );
         THEN( "The triangle normal should be {0,0,1}." ) {
            REQUIRE( triangle.normal[0] == 0.0 );
            REQUIRE( triangle.normal[1] == 0.0 );
            REQUIRE( triangle.normal[2] == 1.0 );
         }
      }
      WHEN( "The normal vector points in negative z-direction and is not normalized." ) {
         constexpr std::array<double, 3> normal = { 0.0, 0.0, -2.34567 };
         StlUtilities::Triangle const triangle( normal, p1, p2, p3 );
         THEN( "The triangle normal should be {0,0,-1}." ) {
            REQUIRE( triangle.normal[0] == 0.0 );
            REQUIRE( triangle.normal[1] == 0.0 );
            REQUIRE( triangle.normal[2] == -1.0 );
         }
      }
      WHEN( "The virtual points are computed" ) {
         constexpr std::array<double, 3> normal = { 0.0, 0.0, 2.34567 };
         StlUtilities::Triangle const triangle( normal, p1, p2, p3 );
         THEN( "The first virtual point should lie at {-1.5, -sin(60), 0.0}." ) {
            REQUIRE( triangle.v[0][0] == Approx( -1.5 ) );
            REQUIRE( triangle.v[0][1] == Approx( -0.8660254038 ) );
            REQUIRE( triangle.v[0][2] == 0.0 );
         }
         THEN( "The second virtual point should lie at {0, 2sin(60), 0.0}." ) {
            REQUIRE( triangle.v[1][0] == Approx( 0.0 ) );
            REQUIRE( triangle.v[1][1] == Approx( 1.732050808 ) );
            REQUIRE( triangle.v[1][2] == 0.0 );
         }
         THEN( "The third virtual point should lie at {1.5, -sin(60), 0.0}." ) {
            REQUIRE( triangle.v[2][0] == Approx( 1.5 ) );
            REQUIRE( triangle.v[2][1] == Approx( -0.8660254038 ) );
            REQUIRE( triangle.v[2][2] == 0.0 );
         }
      }
   }
}

SCENARIO( "The distance to a triangle is computed properly", "[1rank]" ) {
   GIVEN( "A triangle with random non-zero coordinates in all directions and non-normalized edges." ) {
      constexpr std::array<double, 3> p1                   = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> p2                   = { 4.0, 7.0, 8.0 };
      constexpr std::array<double, 3> p3                   = { 7.0, -3.0, -3.5 };
      constexpr std::array<double, 3> normal               = { -7.5, 49.5, -45 };
      constexpr std::array<double, 3> normalized_normal    = { -0.1114141295, 0.7353332544, -0.6684847767 };
      constexpr std::array<double, 3> perpendicular_normal = { 0.9937205673, 0.0752457989, -0.0828100473 };
      constexpr std::array<double, 3> centroid             = { 4.0, 2.0, 2.5 };
      StlUtilities::Triangle const triangle( normal, p1, p2, p3 );
      double distance = std::numeric_limits<double>::max();
      WHEN( "The distance to the triangle is computed for a point along the positive normal in a distance 1.2345678 from the triangle centroid." ) {
         constexpr std::array<double, 3> point = { centroid[0] + normalized_normal[0] * 1.2345678,
                                                   centroid[1] + normalized_normal[1] * 1.2345678,
                                                   centroid[2] + normalized_normal[2] * 1.2345678 };
         THEN( "The distance should be 1.2345678." ) {
            StlUtilities::Voxelization( triangle, point, distance );
            REQUIRE( distance == Approx( 1.2345678 ) );
         }
      }
      WHEN( "The distance to the triangle is computed for a point along the negative normal in a distance 2.537 from the triangle centroid." ) {
         constexpr std::array<double, 3> point = { centroid[0] + normalized_normal[0] * -2.537,
                                                   centroid[1] + normalized_normal[1] * -2.537,
                                                   centroid[2] + normalized_normal[2] * -2.537 };
         THEN( "The distance should be -2.537." ) {
            StlUtilities::Voxelization( triangle, point, distance );
            REQUIRE( distance == Approx( -2.537 ) );
         }
      }
      // This test case indicates the limitations of the algorithm in case the points are too far away from a triangle. If the
      // factor 1.0 is increased to 1.25 there is already a significant deviation. Therefore, the number and size of the triangles
      // must be large enough to ensure that the algorithm works properly.
      WHEN( "The distance to the triangle is computed for a point perpendicular to the positive normal in a distance 2.537." ) {
         constexpr std::array<double, 3> point = { centroid[0] + normalized_normal[0] * 2.537 + perpendicular_normal[0] * 1.0,
                                                   centroid[1] + normalized_normal[1] * 2.537 + perpendicular_normal[1] * 1.0,
                                                   centroid[2] + normalized_normal[2] * 2.537 + perpendicular_normal[2] * 1.0 };
         THEN( "The distance should be 2.537." ) {
            StlUtilities::Voxelization( triangle, point, distance );
            REQUIRE( distance == Approx( 2.537 ) );
         }
      }
      WHEN( "The distance to the triangle is computed for a point distance of 1.0, but a previously assigned distance of 2.5." ) {
         constexpr std::array<double, 3> point = { centroid[0] + normalized_normal[0],
                                                   centroid[1] + normalized_normal[1],
                                                   centroid[2] + normalized_normal[2] };
         distance                              = 2.5;
         THEN( "The distance should be 1.0" ) {
            StlUtilities::Voxelization( triangle, point, distance );
            REQUIRE( distance == Approx( 1.0 ) );
         }
      }
      WHEN( "The distance to the triangle is computed for a point distance of 1.0, but a previously assigned distance of 0.5." ) {
         constexpr std::array<double, 3> point = { centroid[0] + normalized_normal[0],
                                                   centroid[1] + normalized_normal[1],
                                                   centroid[2] + normalized_normal[2] };
         distance                              = 0.5;
         THEN( "The distance should be 0.5" ) {
            StlUtilities::Voxelization( triangle, point, distance );
            REQUIRE( distance == Approx( 0.5 ) );
         }
      }
   }
}
