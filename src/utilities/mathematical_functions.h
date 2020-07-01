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
#ifndef MATHEMATICAL_FUNCTIONS_H
#define MATHEMATICAL_FUNCTIONS_H

#include <algorithm>
#include <cmath>
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Overload of signum function for positive only types ( e.g. unsigned int ).
 */
template<typename T>
typename std::enable_if<std::is_unsigned<T>::value, int>::type
constexpr Signum( T const x ) {
   return T( 0 ) < x;
}

/**
 * @brief Computes the signum, i. e. -1,0,1 of the given input. Not to be mistaken for sgn function which only gives +/- 1.
 */
template<typename T>
typename std::enable_if<std::is_signed<T>::value, int>::type
constexpr Signum( T const x ) {
   return ( T( 0 ) < x ) - ( x < T( 0 ) );
}

/**
 * @brief Computes the sign, i. e. -1,1 of the given input. Note, this functions only gives +/- 1. Zero input results in 1 as result.
 */
template<typename T>
T constexpr Sign( T const x ) {
   return ( T( 0 ) <= x ) - ( x < T( 0 ) );
}

/**
 * @brief Computes the power function for the given input and an integer exponent by repeated multiplication.
 */
template<unsigned int Exponent, typename T>
T constexpr IntegerPow( [[maybe_unused]] T const x ) {
   if constexpr( Exponent == 1 ) {
      return x;
   }
   if constexpr( Exponent > 1 ) {
      return IntegerPow<Exponent-1>(x) * x;
   }
   // else: Exponent == 0
   return T( 1 );
}

/**
 * @brief Computes the floating point consistent sum of three double values. The sum might be not the correctly rounded analytical one,
 *        however, the same result is guaranteed when a, b and c are rearranged.
 */
constexpr double ConsistencyManagedSum( double const a, double const b, double const c ) {
   if constexpr( CC::FUSY() ) {
      double const sum1 = (a + b) + c;
      double const sum2 = (a + c) + b;
      double const sum3 = (b + c) + a;
      return 0.5 * ( std::min({sum1, sum2, sum3}) + std::max({sum1, sum2, sum3}) );
   } else {
      return a + b + c;
   }
}

/**
 * @brief Computes the floating point consistent sum of up to three double values. Dependent on whether the simulation is 1D, 2D or 3D,
 *        only the first, the first two or all three values are considered. For 3D: the sum might be not the correctly rounded analytical one,
 *        however, the same result is guaranteed when a, b and c are rearranged.
 */
constexpr double DimensionAwareConsistencyManagedSum( double const a, double const b, double const c ) {
   if constexpr( CC::DIM() == Dimension::One ) {
      return a;
   }
   if constexpr( CC::DIM() == Dimension::Two ) {
      return a + b;
   }
   return ConsistencyManagedSum( a, b, c );
}

/**
 * @brief Computes the floating point consistent sum of four double values. The sum might be not the correctly rounded analytical one,
 *        however, the same result is guaranteed when a, b, c and d are rearranged.
 */
constexpr double ConsistencyManagedSum( double const a, double const b, double const c, double const d ) {
   if constexpr( CC::FUSY() ) {
      double const sum1 = (a + b) + (c + d);
      double const sum2 = (a + c) + (b + d);
      double const sum3 = (a + d) + (b + c);
      return 0.5 * ( std::min({sum1, sum2, sum3}) + std::max({sum1, sum2, sum3}) );
   } else {
      return a + b + c + d;
   }
}

/**
 * @brief Returns the value given in the array. Auxiliary function.
 */
constexpr double ConsistencyManagedSum( std::array<double, 1> const values ) {
   return values[0];
}

/**
 * @brief Computes the floating point consistent sum of two double values in an array. The sum might be not the correctly rounded analytical one,
 *        however, the same result is guaranteed when a and b are rearranged.
 */
constexpr double ConsistencyManagedSum( std::array<double, 2> const values ) {
   return values[0] + values[1];
}

/**
 * @brief Computes the floating point consistent sum of three double values in an array. The sum might be not the correctly rounded analytical one,
 *        however, the same result is guaranteed when a, b and c are rearranged.
 */
constexpr double ConsistencyManagedSum( std::array<double, 3> const values ) {
   return ConsistencyManagedSum( values[0], values[1], values[2] );
}

/**
 * @brief Computes the floating point consistent sum of up to three double values in an array. Dependent on whether the simulation is 1D, 2D or 3D,
 *        only the first, the first two or all three values are considered. For 3D: the sum might be not the correctly rounded analytical one,
 *        however, the same result is guaranteed when a, b and c are rearranged.
 */
constexpr double DimensionAwareConsistencyManagedSum( std::array<double, 3> const values ) {
   return DimensionAwareConsistencyManagedSum( values[0], values[1], values[2] );
}

/**
 * @brief Computes the floating point consistent sum of four double values in an array. The sum might be not the correctly rounded analytical one,
 *        however, the same result is guaranteed when a, b, c and d are rearranged.
 */
constexpr double ConsistencyManagedSum( std::array<double, 4> const values ) {
   return ConsistencyManagedSum( values[0], values[1], values[2], values[3] );
}

/**
 * @brief Computes the L2 norm of vector 
 * @param vector Vector for which the L2 norm should be computed 
 * @return L2 norm of the vector 
 */
constexpr double L2Norm( std::array<double, 3> const vector ) {
   return std::sqrt( DimensionAwareConsistencyManagedSum( vector[0] * vector[0], vector[1] * vector[1], vector[2] * vector[2] ) );
}

/**
 * @brief Computes the cross product between two vectors 
 * @param vector1, vector2 Two vectors for which the cross product should be computed 
 * @return cross product between the two vectors  
 */
constexpr std::array<double, 3> CrossProduct( std::array<double, 3> const vector1, std::array<double, 3> const vector2 ) {
   return { vector1[1] * vector2[2] - vector1[2] * vector2[1] ,
            vector1[2] * vector2[0] - vector1[0] * vector2[2] ,
            vector1[0] * vector2[1] - vector1[1] * vector2[0] };
}

/**
 * @brief Implementation of the Min-Mod function according to \cite LeVeque1992.
 * @param value_1 The first value.
 * @param value_2 The second value.
 * @return The result of the Min-Mod operation.
 */
inline double MinMod(double const value_1, double const value_2) {
   return 0.5 * ( Signum(value_1) + Signum(value_2) ) * std::min(std::abs(value_1),std::abs(value_2));
}

double GodunovHamiltonian( double const  (&derivatives)[DTI( CC::DIM() )][2], double const  old_levelset_sign );

#endif // MATHEMATICAL_FUNCTIONS_H
