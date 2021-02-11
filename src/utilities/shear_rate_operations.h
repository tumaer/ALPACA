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
#ifndef SHEAR_RATE_UTILITIES_H
#define SHEAR_RATE_UTILITIES_H

#include "stencils/stencil_utilities.h"
#include "utilities/buffer_operations_stencils.h"
#include "utilities/tensor_operations.h"

namespace ShearRateOperations {

   /**
 * @brief Computes the shear rate tensor for a given velocity gradient in the full buffer.
 * @param velocity_gradient Velocity gradient tensor fo which the shear rate is computed.
 * @param shear_rate_tensor Buffer where the shear rate tensor is stored into (indirect return).
 */
   inline void ComputeShearRateTensor( double const ( &velocity_gradient )[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )],
                                       double ( &shear_rate_tensor )[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) {

      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  for( unsigned int c = 0; c < DTI( CC::DIM() ); ++c ) {

                     // Formula: S_ij = 0.5*(du_idx_j + du_jdx_i)
                     shear_rate_tensor[i][j][k][r][c] = 0.5 * ( velocity_gradient[i][j][k][r][c] + velocity_gradient[i][j][k][c][r] );
                  }
               }
            }
         }
      }
   }

   /**
 * @brief Computes the shear rate tensor for a single velocity gradient tensor.
 * @param velocity_gradient Velocity gradient tensor fo which the shear rate is computed.
 * @param shear_rate_tensor tensor where the shear rate tensor is stored into (indirect return).
 */
   inline void ComputeShearRateTensor( std::array<std::array<double, 3>, 3> const& velocity_gradient,
                                       std::array<std::array<double, 3>, 3>& shear_rate_tensor ) {

      for( unsigned int c = 0; c < 3; ++c ) {
         for( unsigned int r = 0; r < 3; ++r ) {

            // Formula: S_ij = 0.5*(du_idx_j + du_jdx_i)
            shear_rate_tensor[c][r] = 0.5 * ( velocity_gradient[c][r] + velocity_gradient[r][c] );
         }
      }
   }

   /**
 * @brief Computes the shear rate for a given velocity vector field.
 * @param u,v,w Velocity buffer for which th shear rate is computed.
 * @param cell_size size of the cell for gradient computation.
 * @param shear_rate Buffer, where the shear rate is stored into (indirect return).
 */
   template<typename DerivativeStencil>
   inline void ComputeShearRate( double const ( &u )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                 double const ( &v )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                 double const ( &w )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                 double const cell_size,
                                 double ( &shear_rate )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {

      /**
    * Description for the positions of the Array:
    * [CC::TCX()]    [CC::TCY()]    [CC::TCZ()]    [DTI(CC::DIM())][DTI(CC::DIM())]
    * Field index x  Field index y  Field index z  Velocity gradient: du_i / dx_j
    */
      double velocity_gradient[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )];

      /**
    * Description for the positions of the Array:
    *  [CC::TCX()]    [CC::TCY()]    [CC::TCZ()]   [DTI(CC::DIM())][DTI(CC::DIM())]
    * Field index x  Field index y  Field index z   Shear-rate: S_ij
    */
      double shear_rate_tensor[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )];

      /**
    * Description for the positions of the Array:
    *  [CC::TCX()]    [CC::TCY()]    [CC::TCZ()]   [DTI(CC::DIM())][DTI(CC::DIM())]
    * Field index x  Field index y  Field index z    Shear-rate^2: S_ij*S_ij
    */
      double shear_rate_tensor_squared[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )];

      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                  for( unsigned int c = 0; c < DTI( CC::DIM() ); ++c ) {
                     velocity_gradient[i][j][k][r][c]         = 0.0;
                     shear_rate_tensor[i][j][k][r][c]         = 0.0;
                     shear_rate_tensor_squared[i][j][k][r][c] = 0.0;
                  }// s
               }   // r
            }      //k
         }         //j
      }            //i

      // Calculate the gradient of the velocity vector
      BO::Stencils::ComputeVectorGradientAtCellCenter<DerivativeStencil>( u, v, w, cell_size, velocity_gradient );
      // compute tensor and squared tensor of the shear_rate
      ComputeShearRateTensor( velocity_gradient, shear_rate_tensor );
      TensorOperations::ComputeSquaredOfTensor( shear_rate_tensor, shear_rate_tensor_squared );

      // Computation of the actual shear rate
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {

               double const second_invariant_of_shear_rate_tensor = TensorOperations::ComputeSecondInvariant( shear_rate_tensor, shear_rate_tensor_squared, i, j, k );

               // Formula: gamma_dot = 2*sqrt(I_2(S_ij))
               shear_rate[i][j][k] = 2.0 * std::sqrt( std::abs( second_invariant_of_shear_rate_tensor ) );
            }
         }
      }
   }

   /**
 * @brief Computes the shear rate for a given velocity vectorial velocity field at a certain position.
 * @param u,v,w Velocity buffer for which th shear rate is computed.
 * @param cell_size size of the cell for gradient computation.
 * @param i, j ,k indices at which position the shear rate should be computed.
 * @return shear rate at the given position.
 */
   template<typename DerivativeStencil>
   inline double ComputeShearRate( double const ( &u )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                   double const ( &v )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                   double const ( &w )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                   double const cell_size,
                                   unsigned int const i, unsigned int const j, unsigned int const k ) {

      // Calculate the velocity gradient
      std::array<std::array<double, 3>, 3> const velocity_gradient   = SU::JacobianMatrix<DerivativeStencil>( u, v, w, i, j, k, cell_size );
      std::array<std::array<double, 3>, 3> shear_rate_tensor         = { { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } } };
      std::array<std::array<double, 3>, 3> shear_rate_tensor_squared = { { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } } };

      // compute tensor and squared tensor of the shear_rate
      ComputeShearRateTensor( velocity_gradient, shear_rate_tensor );
      TensorOperations::ComputeSquaredOfTensor( shear_rate_tensor, shear_rate_tensor_squared );

      // compute the second invariant of the shear rate tensor
      double const second_invariant_of_shear_rate_tensor = TensorOperations::ComputeSecondInvariant( shear_rate_tensor, shear_rate_tensor_squared );

      // retrun the shear_rate
      // Formula: gamma_dot = 2*sqrt(I_2(S_ij))
      return 2.0 * std::sqrt( std::abs( second_invariant_of_shear_rate_tensor ) );
   }

}// namespace ShearRateOperations

#endif// SHEAR_RATE_UTILITIES_H
