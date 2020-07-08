#ifndef TENSOR_OPERATIONS_H
#define TENSOR_OPERATIONS_H

#include "utilities/mathematical_functions.h"
#include "stencils/stencil_utilities.h"

namespace TensorOperations {

   /**
 * @brief Computes the trace of a dim x dim tensor for a full buffer at a certain position
 * @param tensor Buffer containing the dim x dim tensor for each cell 
 * @param i,j,k Index to be used to calcualte the trace 
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
   inline void ComputeSquaredOfTensor( double const ( &tensor )[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )],
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
   inline void ComputeSquaredOfTensor( std::array<std::array<double, 3>, 3> const& tensor,
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

#endif// TENSOR_OPERATIONS_H