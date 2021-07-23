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
#ifndef TENSOR_OPERATIONS_H
#define TENSOR_OPERATIONS_H

#include "utilities/mathematical_functions.h"
#include "stencils/stencil_utilities.h"

namespace TensorOperations {

   /**
 * @brief Computes the trace of a dim x dim tensor for a full buffer at a certain position
 * @param tensor Buffer containing the dim x dim tensor for each cell
 * @param i,j,k Index to be used to calculate the trace
 * @return Trace of tensor at position i,j,k
 */
   constexpr double ComputeTrace( double const ( &tensor )[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )],
                                  unsigned int const i, unsigned int const j, unsigned int const k ) {

      std::array<double, 3> const elements = {
            tensor[i][j][k][0][0],
            CC::DIM() != Dimension::One ? tensor[i][j][k][1][1] : 0.0,
            CC::DIM() == Dimension::Three ? tensor[i][j][k][2][2] : 0.0 };
      return ConsistencyManagedSum( elements );
   }

   /**
 * @brief Computes the second Invariant of a tensor for a full buffer at a given position
 * @param tensor Tensor for which the second invariant should be computed
 * @param tensor_squared squared of the tesnor for which the second invariant should be computed
 * @return Second invarient
 */
   constexpr double ComputeSecondInvariant( double const ( &tensor )[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )],
                                            double const ( &tensor_squared )[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )],
                                            unsigned int const i, unsigned int const j, unsigned int const k ) {

      double const trace_T  = ComputeTrace( tensor, i, j, k );
      double const trace_TT = ComputeTrace( tensor_squared, i, j, k );

      return 0.5 * ( trace_T * trace_T - trace_TT );
   }

   /**
 * @brief Computes the trace of a single dim x dim tensor
 * @param tensor dim x dim tensor for which the trace is computed
 * @return Trace of tensor
 */
   constexpr double ComputeTrace( std::array<std::array<double, 3>, 3> const& tensor ) {

      return ConsistencyManagedSum( tensor[0][0], tensor[1][1], tensor[2][2] );
   }

   /**
 * @brief Computes the second Invariant of a single tensor
 * @param tensor Tensor for which the second invariant should be computed
 * @param tensor_squared Squared of the tensor for which the second invariant should be computed
 * @return Second invarient of tensor
 */
   constexpr double ComputeSecondInvariant( std::array<std::array<double, 3>, 3> const& tensor, std::array<std::array<double, 3>, 3> const& tensor_squared ) {

      double const trace_T  = ComputeTrace( tensor );
      double const trace_TT = ComputeTrace( tensor_squared );

      return 0.5 * ( trace_T * trace_T - trace_TT );
   }

   /**
 * @brief Computes the squared of a tensor for a full buffer
 * @param tensor Tensor for which the squared should be calculated
 * @param tensor_squared Squared of the tensor (as indirect return)
 */
   constexpr void ComputeSquaredOfTensor( double const ( &tensor )[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )],
                                          double ( &tensor_squared )[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) {

      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  for( unsigned int c = 0; c < DTI( CC::DIM() ); ++c ) {

                     std::array<double, 3> const elements = {
                           tensor[i][j][k][r][0] * tensor[i][j][k][0][c],
                           CC::DIM() != Dimension::One ? tensor[i][j][k][r][1] * tensor[i][j][k][1][c] : 0.0,
                           CC::DIM() == Dimension::Three ? tensor[i][j][k][r][2] * tensor[i][j][k][2][c] : 0.0 };

                     tensor_squared[i][j][k][r][c] = ConsistencyManagedSum( elements );
                  }
               }
            }
         }
      }
   }

   /**
 * @brief Computes the squared of a tensor for a single tensor
 * @param tensor tensor for which the squared is computed
 * @param tensor_squared tensor where the squared is stored into (inidrect return)
 */
   constexpr void ComputeSquaredOfTensor( std::array<std::array<double, 3>, 3> const& tensor,
                                          std::array<std::array<double, 3>, 3>& tensor_squared ) {

      for( unsigned int r = 0; r < 3; ++r ) {
         for( unsigned int c = 0; c < 3; ++c ) {
            std::array<double, 3> const elements = {
                  tensor[0][r] * tensor[c][0], tensor[1][r] * tensor[c][1], tensor[2][r] * tensor[c][2] };
            tensor_squared[c][r] = ConsistencyManagedSum( elements );
         }
      }
   }
}// namespace TensorOperations

namespace TO = TensorOperations;

#endif// TENSOR_OPERATIONS_H
