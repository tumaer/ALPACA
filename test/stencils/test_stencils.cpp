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
#include <algorithm>

#include "stencils/stencil_utilities.h"

template<typename S>
void TestStencilSizes( S const& stencil, unsigned int const stencil_size, unsigned int const downstream_size ) {
   REQUIRE( stencil.StencilSize() == stencil_size );
   REQUIRE( stencil.DownstreamStencilSize() == downstream_size );
}

template<typename S>
void TestStencilType( S const& stencil, StencilType const type ) {
   REQUIRE( stencil.GetStencilType() == type );
}

template<typename S>
void TestStencilParameters( S const& stencil, StencilType const type, unsigned int const stencil_size, unsigned int const downstream_size ) {
   THEN( "Stencil size, upstream stencil size and stencil type are correct" ) {
      TestStencilType( stencil, type );
      TestStencilSizes( stencil, stencil_size, downstream_size );
   }
}

template<typename S>
void TestUnityOnUnitArrayWithMargin( double const allowed_margin ) {
   WHEN( "Stencil is applied to unit-array (with cell size one)" ) {
      std::array<double, S::StencilSize()> unit_array;
      unit_array.fill( 1.0 );
      constexpr double cell_size = 1.0;
      THEN( "The resulting value is unity" ) {
         REQUIRE( StencilUtilities::Reconstruction<S,StencilProperty::UpwindLeft,double>( unit_array, cell_size ) == Approx( 1.0 ).margin( allowed_margin ) );
         REQUIRE( StencilUtilities::Reconstruction<S,StencilProperty::UpwindRight,double>( unit_array, cell_size ) == Approx( 1.0 ).margin( allowed_margin ) );
      }
   }
}

template<typename S>
void TestWenoValuesLeftAndRightOfStepWithMargin( double const allowed_margin ) {
   WHEN( "Stencil is applied to unit step-function array (with cell size one)" ) {
      unsigned int const size = S::StencilSize();
      std::array<double, size> step_array;
      step_array.fill( 1.0 );
      std::fill_n( step_array.begin(), size / 2 , 0.0 );
      constexpr double cell_size = 1.0;
      THEN( "The left or right value of the step is returned - according to upwinding direction" ) {
         REQUIRE( StencilUtilities::Reconstruction<S,StencilProperty::UpwindLeft,double>( step_array, cell_size ) == Approx( 0.0 ).margin( allowed_margin ) );
         REQUIRE( StencilUtilities::Reconstruction<S,StencilProperty::UpwindRight,double>( step_array, cell_size ) == Approx( 1.0 ).margin( allowed_margin ) );
      }
   }
}

template<typename S>
void TestZeronessOnUnitArrayWithMargin( double const allowed_margin ) {
   WHEN( "Stencil is applied to unit-array (with cell size one)" ) {
      unsigned int const size = S::StencilSize();
      std::array<double, size> unit_array;
      unit_array.fill( 1.0 );
      constexpr double cell_size = 1.0;
      THEN( "The resulting value is zero" ) {
         REQUIRE( SU::Derivative<S,StencilProperty::UpwindLeft>( unit_array, cell_size ) == Approx( 0.0 ).margin( allowed_margin ) );
         REQUIRE( SU::Derivative<S,StencilProperty::UpwindRight>( unit_array, cell_size ) == Approx( 0.0 ).margin( allowed_margin ) );
      }
   }
}

SCENARIO( "Reconstruction stencil correctness", "[1rank]" ) {
   GIVEN( "A WENO-3 reconstruction stencil" ) {
      constexpr auto weno3 = WENO3();
      TestStencilParameters( weno3, StencilType::Reconstruction, 4, 1 );
      TestUnityOnUnitArrayWithMargin<WENO3>( 1e-16 );
      TestWenoValuesLeftAndRightOfStepWithMargin<WENO3>( 1e-16 );
   }
   GIVEN( "A WENO-5 reconstruction stencil" ) {
      constexpr auto weno5 = WENO5();
      TestStencilParameters( weno5, StencilType::Reconstruction, 6, 2 );
      TestUnityOnUnitArrayWithMargin<WENO5>( 1e-16 );
      TestWenoValuesLeftAndRightOfStepWithMargin<WENO5>( 1e-11 );
   }
   GIVEN( "A WENO-5-Z reconstruction stencil" ) {
      constexpr auto weno5z = WENO5Z();
      TestStencilParameters( weno5z, StencilType::Reconstruction, 6, 2 );
      TestUnityOnUnitArrayWithMargin<WENO5Z>( 1e-16 );
      TestWenoValuesLeftAndRightOfStepWithMargin<WENO5Z>( 1e-15 );
   }
   GIVEN( "A WENO-7 reconstruction stencil" ) {
      constexpr auto weno7 = WENO7();
      TestStencilParameters( weno7, StencilType::Reconstruction, 8, 3 );
      TestUnityOnUnitArrayWithMargin<WENO7>( 1e-16 );
      TestWenoValuesLeftAndRightOfStepWithMargin<WENO7>( 1e-16 );
   }
   GIVEN( "A WENO-9 reconstruction stencil" ) {
      constexpr auto weno9 = WENO9();
      TestStencilParameters( weno9, StencilType::Reconstruction, 10, 4 );
      TestUnityOnUnitArrayWithMargin<WENO9>( 1e-16 );
      TestWenoValuesLeftAndRightOfStepWithMargin<WENO9>( 1e-16 );
   }

   GIVEN( "A WENO-AO53 reconstruction stencil" ) {
      constexpr auto wenoao53 = WENOAO53();
      TestStencilParameters( wenoao53, StencilType::Reconstruction, 6, 2 );
      TestUnityOnUnitArrayWithMargin<WENOAO53>( 1e-16 );
      TestWenoValuesLeftAndRightOfStepWithMargin<WENOAO53>( 1e-16 );
   }

   GIVEN( "A WENO-5HM reconstruction stencil" ) {
      constexpr auto weno5hm = WENO5HM();
      TestStencilParameters( weno5hm, StencilType::Reconstruction, 6, 2 );
      TestUnityOnUnitArrayWithMargin<WENO5HM>( 1e-16 );
      TestWenoValuesLeftAndRightOfStepWithMargin<WENO5HM>( 1e-15 );
   }

   GIVEN( "A WENO-CU6 reconstruction stencil" ) {
      constexpr auto wenocu6 = WENOCU6();
      TestStencilParameters( wenocu6, StencilType::Reconstruction, 6, 2 );
      TestUnityOnUnitArrayWithMargin<WENOCU6>( 1e-16 );
      TestWenoValuesLeftAndRightOfStepWithMargin<WENOCU6>( 1e-8 );
   }

   GIVEN( "A TENO-5 reconstruction stencil" ) {
      constexpr auto teno5 = TENO5();
      TestStencilParameters( teno5, StencilType::Reconstruction, 6, 2 );
      TestUnityOnUnitArrayWithMargin<TENO5>( 1e-16 );
      TestWenoValuesLeftAndRightOfStepWithMargin<TENO5>( 1e-16 );
   }

   GIVEN( "A first order reconstruction stencil" ) {
      constexpr auto first_order = FirstOrder();
      TestStencilParameters( first_order, StencilType::Reconstruction, 2, 0 );
      TestUnityOnUnitArrayWithMargin<FirstOrder>( 0.0 );
      //First order must give WENO-like upwind value ...
      TestWenoValuesLeftAndRightOfStepWithMargin<FirstOrder>( 0.0 );
   }

   GIVEN( "A fourth-order central reconstruction stencil" ) {
      constexpr auto fourth_order = FourthOrderCentral();
      TestStencilParameters( fourth_order, StencilType::Reconstruction, 4, 1 );
      TestUnityOnUnitArrayWithMargin<FourthOrderCentral>( 1e-16 );
      //Central (as the name suggests) does not follow WENO-like upwinding
      WHEN( "Stencil is applied to unit step-function array (with cell size one)" ) {
         unsigned int const size = FourthOrderCentral::StencilSize();
         std::array<double, size> step_array;
         step_array.fill( 1.0 );
         std::fill_n( step_array.begin(), size / 2 , 0.0 );
         constexpr double cell_size = 1.0;
         THEN( "0.5 is returned regardless the upwinding direction" ) {
            REQUIRE( StencilUtilities::Reconstruction<FourthOrderCentral,StencilProperty::UpwindLeft, double>( step_array, cell_size ) == Approx( 0.5 ).margin( 1e-16 ) );
            REQUIRE( StencilUtilities::Reconstruction<FourthOrderCentral,StencilProperty::UpwindRight, double>( step_array, cell_size ) == Approx( 0.5 ).margin( 1e-16 ) );
         }
      }
   }
}

SCENARIO( "Derivative stencil correctness", "[1rank]" ) {
   GIVEN( "A HOUC-5 derivative stencil" ) {
      constexpr auto houc5 = HOUC5();
      TestStencilParameters( houc5, StencilType::Derivative, 7, 3 );
      TestZeronessOnUnitArrayWithMargin<HOUC5>( 1e-16 );
      WHEN( "Stencil is applied to unit step-function array (with cell size one)" ) {
         unsigned int const size = houc5.StencilSize();
         std::array<double, size> step_array;
         step_array.fill( 1.0 );
         std::fill_n( step_array.begin(), size / 2 , 0.0 );
         constexpr double cell_size = 1.0;
         THEN( "0.78333 is returned for upwind right direction and 0.45 for upwind left" ) {
            REQUIRE( SU::Derivative<HOUC5,StencilProperty::UpwindLeft>( step_array, cell_size ) == Approx( 0.7833333333333333333333333333333 ).margin( 1e-16 ) );
            REQUIRE( SU::Derivative<HOUC5,StencilProperty::UpwindRight>( step_array, cell_size ) == Approx( 0.45 ).margin( 1e-16 ) );
         }
      }   
   }

   GIVEN( "A fourth order central difference derivative stencil" ) {
      constexpr auto fourth_order_central_difference = FourthOrderCentralDifference();
      TestStencilParameters( fourth_order_central_difference, StencilType::Derivative, 5, 2 );
      TestZeronessOnUnitArrayWithMargin<FourthOrderCentralDifference>( 1e-16 );
      WHEN( "Stencil is applied to unit step-function array (with cell size one)" ) {
         unsigned int const size = fourth_order_central_difference.StencilSize();
         std::array<double, size> step_array;
         step_array.fill( 1.0 );
         std::fill_n( step_array.begin(), size / 2 , 0.0 );
         constexpr double cell_size = 1.0;
         THEN( "0.58333 is returned for the central difference direction" ) {
            REQUIRE( SU::Derivative<FourthOrderCentralDifference,StencilProperty::Central>( step_array, cell_size ) == Approx( 0.583333333333333333333 ).margin( 1e-16 ) );
         }
      }   
   }

   GIVEN( "A fourth order cell face derivative stencil" ) {
      constexpr auto fourth_order_cell_face = FourthOrderCellFace();
      TestStencilParameters( fourth_order_cell_face, StencilType::Derivative, 4, 1 );
      TestZeronessOnUnitArrayWithMargin<FourthOrderCellFace>( 1e-16 );
      WHEN( "Stencil is applied to unit step-function array (with cell size one)" ) {
         unsigned int const size = fourth_order_cell_face.StencilSize();
         std::array<double, size> step_array;
         step_array.fill( 1.0 );
         std::fill_n( step_array.begin(), size / 2 , 0.0 );
         constexpr double cell_size = 1.0;
         THEN( "1.0833333 is returned for the upwind left direction" ) {
            REQUIRE( SU::Derivative<FourthOrderCellFace,StencilProperty::UpwindLeft>( step_array, cell_size ) == Approx( 1.083333333333333333333 ).margin( 1e-16 ) );
            REQUIRE( SU::Derivative<FourthOrderCellFace,StencilProperty::UpwindRight>( step_array, cell_size ) == Approx( 1.083333333333333333333 ).margin( 1e-16 ) );
         }
      }   
   }

   GIVEN( "A central difference derivative stencil" ) {
      constexpr auto central_difference = CentralDifference();
      TestStencilParameters( central_difference, StencilType::Derivative, 3, 1 );
      TestZeronessOnUnitArrayWithMargin<CentralDifference>( 1e-16 );
      WHEN( "Stencil is applied to unit step-function array (with cell size one)" ) {
         unsigned int const size = central_difference.StencilSize();
         std::array<double, size> step_array;
         step_array.fill( 1.0 );
         std::fill_n( step_array.begin(), size / 2 , 0.0 );
         constexpr double cell_size = 1.0;
         THEN( "0.5 is returned for the central difference direction" ) {
            REQUIRE( SU::Derivative<CentralDifference,StencilProperty::Central>( step_array, cell_size ) == Approx( 0.5 ).margin( 1e-16 ) );
         }
      }   
   }
}