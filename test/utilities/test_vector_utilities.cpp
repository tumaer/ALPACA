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

#include "utilities/vector_utilities.h"

SCENARIO( "The scalar multiplication works properly", "[1rank]" ) {
   GIVEN( "A vector with three distinct variables." ) {
      constexpr std::array<double, 3> vector = { 1.0, 2.0, 3.0 };
      WHEN( "The multiplication factor is 0.0." ) {
         std::array<double, 3> const result = VU::Multiply( vector, 0.0 );
         THEN( "The resulting vector should be {0.0, 0.0, 0.0}" ) {
            for( double const& res : result ) {
               REQUIRE( res == 0.0 );
            }
         }
      }
      WHEN( "The multiplication factor is 1.5." ) {
         std::array<double, 3> const result = VU::Multiply( vector, 1.5 );
         THEN( "The resulting vector should be {1.5, 3.0, 4.5}" ) {
            REQUIRE( result[0] == 1.5 );
            REQUIRE( result[1] == 3.0 );
            REQUIRE( result[2] == 4.5 );
         }
      }
      WHEN( "The multiplication factor is -1.25." ) {
         std::array<double, 3> const result = VU::Multiply( vector, -1.25 );
         THEN( "The resulting vector should be {-1.25, -2.5, -3.75}" ) {
            REQUIRE( result[0] == -1.25 );
            REQUIRE( result[1] == -2.5 );
            REQUIRE( result[2] == -3.75 );
         }
      }
   }
}

SCENARIO( "The L2 norm works properly", "[1rank]" ) {
   GIVEN( "A zero vector." ) {
      constexpr std::array<double, 3> vector = { 0.0, 0.0, 0.0 };
      THEN( "The L2 norm should be zero" ) {
         REQUIRE( VU::L2Norm( vector ) == 0.0 );
      }
   }
   GIVEN( "A vector with unity in the x component." ) {
      constexpr std::array<double, 3> vector = { 1.0, 0.0, 0.0 };
      THEN( "The L2 norm should be unity" ) {
         REQUIRE( VU::L2Norm( vector ) == 1.0 );
      }
   }
   GIVEN( "A vector with unity in the y component." ) {
      constexpr std::array<double, 3> vector = { 0.0, 1.0, 0.0 };
      THEN( "The L2 norm should be unity" ) {
         REQUIRE( VU::L2Norm( vector ) == 1.0 );
      }
   }
   GIVEN( "A vector with unity in the z component." ) {
      constexpr std::array<double, 3> vector = { 0.0, 0.0, 1.0 };
      THEN( "The L2 norm should be unity" ) {
         REQUIRE( VU::L2Norm( vector ) == 1.0 );
      }
   }
   GIVEN( "A vector with three distinct components." ) {
      constexpr std::array<double, 3> vector = { 1.25, 1.3, -0.75 };
      THEN( "The L2 norm should be 1.953202498" ) {
         REQUIRE( VU::L2Norm( vector ) == Approx( 1.9532024985 ) );
      }
   }
}

SCENARIO( "The cross product works properly", "[1rank]" ) {
   GIVEN( "A non-zero v1={1,2,3} and zero vector v2={0,0,0}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> vector2 = { 0.0, 0.0, 0.0 };
      WHEN( "The zero vector is the first input argument. " ) {
         constexpr std::array<double, 3> result = VU::CrossProduct( vector1, vector2 );
         THEN( "The cross product should be zero" ) {
            for( double const& res : result ) {
               REQUIRE( res == 0.0 );
            }
         }
      }
      WHEN( "The zero vector is the second input argument. " ) {
         constexpr std::array<double, 3> result = VU::CrossProduct( vector2, vector1 );
         THEN( "The cross product should be zero" ) {
            for( double const& res : result ) {
               REQUIRE( res == 0.0 );
            }
         }
      }
   }
   GIVEN( "The unit vectors e_x = {1,0,0} and e_y={0,1,0}." ) {
      constexpr std::array<double, 3> e_x = { 1.0, 0.0, 0.0 };
      constexpr std::array<double, 3> e_y = { 0.0, 1.0, 0.0 };
      constexpr std::array<double, 3> e_z = VU::CrossProduct( e_x, e_y );
      THEN( "The cross product should be e_z = {0,0,1}" ) {
         REQUIRE( e_z[0] == 0.0 );
         REQUIRE( e_z[1] == 0.0 );
         REQUIRE( e_z[2] == 1.0 );
      }
   }
   GIVEN( "Two non-zero vectors v1 = {1,2,3} and v2={4,5,6}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> vector2 = { 4.0, 5.0, 6.0 };
      WHEN( "The cross product is computed." ) {
         constexpr std::array<double, 3> result = VU::CrossProduct( vector1, vector2 );
         THEN( "The resulting vector should be {-3,6-3}" ) {
            REQUIRE( result[0] == -3.0 );
            REQUIRE( result[1] == 6.0 );
            REQUIRE( result[2] == -3.0 );
         }
      }
      WHEN( "The vectors are permuted as input arguments." ) {
         constexpr std::array<double, 3> result1 = VU::CrossProduct( vector1, vector2 );
         constexpr std::array<double, 3> result2 = VU::CrossProduct( vector2, vector1 );
         THEN( "The cross product should be the negative of the original." ) {
            for( unsigned int i = 0; i < 3; i++ ) {
               REQUIRE( result1[i] == -result2[i] );
            }
         }
      }
   }
}

SCENARIO( "The dot product works properly", "[1rank]" ) {
   GIVEN( "A non-zero v1={1,2,3} and zero vector v2={0,0,0}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> vector2 = { 0.0, 0.0, 0.0 };
      WHEN( "The zero vector is the first input argument. " ) {
         THEN( "The dot product should be zero" ) {
            REQUIRE( VU::DotProduct( vector1, vector2 ) == 0.0 );
         }
      }
      WHEN( "The zero vector is the second input argument. " ) {
         THEN( "The dot product should be zero" ) {
            REQUIRE( VU::DotProduct( vector1, vector2 ) == 0.0 );
         }
      }
   }
   GIVEN( "The unit vectors e_x = {1,0,0} and e_y={0,1,0}." ) {
      constexpr std::array<double, 3> e_x = { 1.0, 0.0, 0.0 };
      constexpr std::array<double, 3> e_y = { 0.0, 1.0, 0.0 };
      THEN( "The dot product should be zero" ) {
         REQUIRE( VU::DotProduct( e_x, e_y ) == 0.0 );
      }
   }
   GIVEN( "Two non-zero vectors v1 = {1,2,3} and v2={4,5,-6}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> vector2 = { 4.0, 5.0, -6.0 };
      WHEN( "The dot product is computed." ) {
         THEN( "The dot product should be 32" ) {
            REQUIRE( VU::DotProduct( vector1, vector2 ) == -4.0 );
         }
      }
      WHEN( "The vectors are permuted as input arguments." ) {
         THEN( "The dot products should be the same." ) {
            REQUIRE( VU::DotProduct( vector1, vector2 ) == VU::DotProduct( vector2, vector1 ) );
         }
      }
   }
}

SCENARIO( "The normalization works properly", "[1rank]" ) {
   GIVEN( "The unit vector = {1,0,0}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 0.0, 0.0 };
      THEN( "The resulting vector should be {1,0,0}" ) {
         std::array<double, 3> const result = VU::Normalize( vector1 );
         REQUIRE( result[0] == 1.0 );
         REQUIRE( result[1] == 0.0 );
         REQUIRE( result[2] == 0.0 );
      }
   }
   GIVEN( "The vector = {1,1,1}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 1.0, 1.0 };
      std::array<double, 3> const result      = VU::Normalize( vector1 );
      THEN( "Each component of the resulting vector should be 0.5773502692" ) {
         for( unsigned int i = 0; i < 3; i++ ) {
            REQUIRE( result[i] == Approx( 0.5773502692 ) );
         }
      }
      THEN( "The length of the vector should be one" ) {
         REQUIRE( std::sqrt( result[0] * result[0] + result[1] * result[1] + result[2] * result[2] ) == 1.0 );
      }
   }
}

SCENARIO( "The difference works properly", "[1rank]" ) {
   GIVEN( "A single non-zero vector v1 = {1,2,3}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 2.0, 3.0 };
      std::array<double, 3> const result      = VU::Difference( vector1, vector1 );
      THEN( "The difference of the vector and itself should be the zero vector." ) {
         for( double const& res : result ) {
            REQUIRE( res == 0.0 );
         }
      }
   }
   GIVEN( "A non-zero v1={1,2,3} and zero vector v2={0,0,0}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> vector2 = { 0.0, 0.0, 0.0 };
      WHEN( "The zero vector is the first input argument. " ) {
         std::array<double, 3> const result = VU::Difference( vector2, vector1 );
         THEN( "The diffence should be {1,2,3}" ) {
            REQUIRE( result[0] == 1.0 );
            REQUIRE( result[1] == 2.0 );
            REQUIRE( result[2] == 3.0 );
         }
      }
      WHEN( "The zero vector is the second input argument it should be the negative than zero as first argument. " ) {
         std::array<double, 3> const result1 = VU::Difference( vector2, vector1 );
         std::array<double, 3> const result2 = VU::Difference( vector1, vector2 );
         THEN( "The difference should be the negative of the original." ) {
            for( unsigned int i = 0; i < 3; i++ ) {
               REQUIRE( result1[i] == -result2[i] );
            }
         }
      }
   }
   GIVEN( "Two non-zero vectors v1 = {1,2,3} and v2={4,5,6}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> vector2 = { 4.0, 5.0, 6.0 };
      std::array<double, 3> const result      = VU::Difference( vector1, vector2 );
      THEN( "The difference should be 3 in all components" ) {
         for( double const& res : result ) {
            REQUIRE( res == 3.0 );
         }
      }
   }
}

SCENARIO( "The distance works properly", "[1rank]" ) {
   GIVEN( "A single non-zero vector v1 = {1,2,3}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 2.0, 3.0 };
      THEN( "The distance of the vector and itself should be zero" ) {
         REQUIRE( VU::Distance( vector1, vector1 ) == 0.0 );
      }
   }
   GIVEN( "Two non-zero vectors v1 = {1,2,3} and v2={4,5,6}." ) {
      constexpr std::array<double, 3> vector1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> vector2 = { 4.0, 5.0, 6.0 };
      WHEN( "The distance is computed with v1 as first input argument." ) {
         THEN( "The distance should be 5.196152423" ) {
            REQUIRE( VU::Distance( vector1, vector2 ) == Approx( 5.196152423 ) );
         }
      }
      WHEN( "The input arguments are permuted." ) {
         THEN( "The distance should be the same" ) {
            REQUIRE( VU::Distance( vector1, vector2 ) == VU::Distance( vector2, vector1 ) );
         }
      }
   }
}

SCENARIO( "The scalar triple product works properly", "[1rank]" ) {
   GIVEN( "The unit vectors e_x = {1,0,0}, e_y={0,1,0} and e_z={0,0,1}." ) {
      constexpr std::array<double, 3> e_x = { 1.0, 0.0, 0.0 };
      constexpr std::array<double, 3> e_y = { 0.0, 1.0, 0.0 };
      constexpr std::array<double, 3> e_z = { 0.0, 0.0, 1.0 };
      THEN( "The triple product should be 1.0" ) {
         REQUIRE( VU::ScalarTripleProduct( e_x, e_y, e_z ) == 1.0 );
      }
   }
   GIVEN( "Two arbitrary vectors v1 = {1,2,3}, v2={4,5,6} and a zero vector={0,0,0}." ) {
      constexpr std::array<double, 3> v1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> v2 = { 4.0, 5.0, 6.0 };
      constexpr std::array<double, 3> v3 = { 0.0, 0.0, 0.0 };
      THEN( "The triple product should be 0.0 for all permutations" ) {
         REQUIRE( VU::ScalarTripleProduct( v1, v2, v3 ) == 0.0 );
         REQUIRE( VU::ScalarTripleProduct( v1, v3, v2 ) == 0.0 );
         REQUIRE( VU::ScalarTripleProduct( v2, v1, v3 ) == 0.0 );
         REQUIRE( VU::ScalarTripleProduct( v2, v3, v1 ) == 0.0 );
         REQUIRE( VU::ScalarTripleProduct( v3, v1, v2 ) == 0.0 );
         REQUIRE( VU::ScalarTripleProduct( v3, v2, v1 ) == 0.0 );
      }
   }
   GIVEN( "Three non-zero unit vectors v1 = {1,2,3}, v2={4,5,6} and v3={7,8,9}." ) {
      constexpr std::array<double, 3> v1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> v2 = { 4.0, 5.0, 6.0 };
      constexpr std::array<double, 3> v3 = { 7.0, 8.0, 9.0 };
      THEN( "The triple product should be 0.0" ) {
         REQUIRE( VU::ScalarTripleProduct( v1, v2, v3 ) == 0.0 );
      }
   }
   GIVEN( "Three non-zero unit vectors v1 = {1,2,3}, v2={4,5,6} and v3={7,8,10}." ) {
      constexpr std::array<double, 3> v1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> v2 = { 4.0, 5.0, 6.0 };
      constexpr std::array<double, 3> v3 = { 7.0, 8.0, 10.0 };
      WHEN( "The triple product is computed." ) {
         THEN( "The triple product should be -3.0" ) {
            REQUIRE( VU::ScalarTripleProduct( v1, v2, v3 ) == -3.0 );
         }
      }
      WHEN( "The three vectors undergo cyclic permutation." ) {
         THEN( "The triple product should be the same for all" ) {
            REQUIRE( VU::ScalarTripleProduct( v1, v2, v3 ) == VU::ScalarTripleProduct( v3, v1, v2 ) );
            REQUIRE( VU::ScalarTripleProduct( v1, v2, v3 ) == VU::ScalarTripleProduct( v2, v3, v1 ) );
         }
      }
      WHEN( "The two vectors are permuted." ) {
         THEN( "The triple product should be the negative of the original." ) {
            REQUIRE( VU::ScalarTripleProduct( v1, v2, v3 ) == -VU::ScalarTripleProduct( v2, v1, v3 ) );
            REQUIRE( VU::ScalarTripleProduct( v1, v2, v3 ) == -VU::ScalarTripleProduct( v3, v2, v1 ) );
            REQUIRE( VU::ScalarTripleProduct( v1, v2, v3 ) == -VU::ScalarTripleProduct( v1, v3, v2 ) );
         }
      }
   }
}

SCENARIO( "The vector triple product works properly", "[1rank]" ) {
   GIVEN( "The unit vectors e_x = {1,0,0}, e_y={0,1,0} and e_z={0,0,1}." ) {
      constexpr std::array<double, 3> e_x    = { 1.0, 0.0, 0.0 };
      constexpr std::array<double, 3> e_y    = { 0.0, 1.0, 0.0 };
      constexpr std::array<double, 3> e_z    = { 0.0, 0.0, 1.0 };
      constexpr std::array<double, 3> result = VU::VectorTripleProduct( e_x, e_y, e_z );
      THEN( "The triple product should be 0.0 in all components" ) {
         for( double const& res : result ) {
            REQUIRE( res == 0.0 );
         }
      }
   }
   GIVEN( "Three non-zero unit vectors v1 = {1,2,3}, v2={4,5,6} and v3={7,8,9}." ) {
      constexpr std::array<double, 3> v1 = { 1.0, 2.0, 3.0 };
      constexpr std::array<double, 3> v2 = { 4.0, 5.0, 6.0 };
      constexpr std::array<double, 3> v3 = { 7.0, 8.0, 9.0 };
      WHEN( "The triple product is computed." ) {
         constexpr std::array<double, 3> result = VU::VectorTripleProduct( v1, v2, v3 );
         THEN( "The triple product should be {78, 6, -66}" ) {
            REQUIRE( result[0] == 78.0 );
            REQUIRE( result[1] == 6.0 );
            REQUIRE( result[2] == -66.0 );
         }
      }
      WHEN( "The first two vectors are permuted." ) {
         constexpr std::array<double, 3> result1 = VU::VectorTripleProduct( v1, v2, v3 );
         constexpr std::array<double, 3> result2 = VU::VectorTripleProduct( v2, v1, v3 );
         THEN( "The triple product should be the negative of the original." ) {
            for( unsigned int i = 0; i < 3; i++ ) {
               REQUIRE( result1[i] == -result2[i] );
            }
         }
      }
   }
}
