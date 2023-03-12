//===-------------------- shear_rate_operations.h -------------------------===//
//
//                                 ALPACA
//
// Part of ALPACA, under the GNU General Public License as published by
// the Free Software Foundation version 3.
// SPDX-License-Identifier: GPL-3.0-only
//
// If using this code in an academic setting, please cite the following:
// @article{hoppe2022parallel,
//  title={A parallel modular computing environment for three-dimensional
//  multiresolution simulations of compressible flows},
//  author={Hoppe, Nils and Adami, Stefan and Adams, Nikolaus A},
//  journal={Computer Methods in Applied Mechanics and Engineering},
//  volume={391},
//  pages={114486},
//  year={2022},
//  publisher={Elsevier}
// }
//
//===----------------------------------------------------------------------===//
#ifndef SHEAR_RATE_UTILITIES_H
#define SHEAR_RATE_UTILITIES_H

#include "stencils/stencil_utilities.h"
#include "utilities/buffer_operations_stencils.h"
#include "utilities/tensor_operations.h"

namespace ShearRateOperations {

/**
 * @brief Computes the shear rate tensor for a given velocity gradient in the
 * full buffer.
 * @param velocity_gradient Velocity gradient tensor fo which the shear rate is
 * computed.
 * @param shear_rate_tensor Buffer where the shear rate tensor is stored into
 * (indirect return).
 */
inline void ComputeShearRateTensor(
    double const (&velocity_gradient)[CC::TCX()][CC::TCY()][CC::TCZ()]
                                     [DTI(CC::DIM())][DTI(CC::DIM())],
    double (&shear_rate_tensor)[CC::TCX()][CC::TCY()][CC::TCZ()][DTI(CC::DIM())]
                               [DTI(CC::DIM())]) {

  for (unsigned int i = 0; i < CC::TCX(); ++i) {
    for (unsigned int j = 0; j < CC::TCY(); ++j) {
      for (unsigned int k = 0; k < CC::TCZ(); ++k) {
        for (unsigned int r = 0; r < DTI(CC::DIM()); ++r) {
          for (unsigned int c = 0; c < DTI(CC::DIM()); ++c) {

            // Formula: S_ij = 0.5*(du_idx_j + du_jdx_i)
            shear_rate_tensor[i][j][k][r][c] =
                0.5 * (velocity_gradient[i][j][k][r][c] +
                       velocity_gradient[i][j][k][c][r]);
          }
        }
      }
    }
  }
}

/**
 * @brief Computes the shear rate tensor for a single velocity gradient tensor.
 * @param velocity_gradient Velocity gradient tensor fo which the shear rate is
 * computed.
 * @param shear_rate_tensor tensor where the shear rate tensor is stored into
 * (indirect return).
 */
inline void ComputeShearRateTensor(
    std::array<std::array<double, 3>, 3> const &velocity_gradient,
    std::array<std::array<double, 3>, 3> &shear_rate_tensor) {

  for (unsigned int c = 0; c < 3; ++c) {
    for (unsigned int r = 0; r < 3; ++r) {

      // Formula: S_ij = 0.5*(du_idx_j + du_jdx_i)
      shear_rate_tensor[c][r] =
          0.5 * (velocity_gradient[c][r] + velocity_gradient[r][c]);
    }
  }
}

/**
 * @brief Computes the shear rate for a given velocity vector field.
 * @param u,v,w Velocity buffer for which th shear rate is computed.
 * @param cell_size size of the cell for gradient computation.
 * @param shear_rate Buffer, where the shear rate is stored into (indirect
 * return).
 */
template <typename DerivativeStencil>
inline void
ComputeShearRate(double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()],
                 double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()],
                 double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()],
                 double const cell_size,
                 double (&shear_rate)[CC::TCX()][CC::TCY()][CC::TCZ()]) {

  /**
   * Description for the positions of the Array:
   * [CC::TCX()]    [CC::TCY()]    [CC::TCZ()] [DTI(CC::DIM())][DTI(CC::DIM())]
   * Field index x  Field index y  Field index z  Velocity gradient: du_i / dx_j
   */
  double velocity_gradient[CC::TCX()][CC::TCY()][CC::TCZ()][DTI(CC::DIM())]
                          [DTI(CC::DIM())];

  /**
   * Description for the positions of the Array:
   *  [CC::TCX()]    [CC::TCY()]    [CC::TCZ()] [DTI(CC::DIM())][DTI(CC::DIM())]
   * Field index x  Field index y  Field index z   Shear-rate: S_ij
   */
  double shear_rate_tensor[CC::TCX()][CC::TCY()][CC::TCZ()][DTI(CC::DIM())]
                          [DTI(CC::DIM())];

  /**
   * Description for the positions of the Array:
   *  [CC::TCX()]    [CC::TCY()]    [CC::TCZ()] [DTI(CC::DIM())][DTI(CC::DIM())]
   * Field index x  Field index y  Field index z    Shear-rate^2: S_ij*S_ij
   */
  double shear_rate_tensor_squared[CC::TCX()][CC::TCY()][CC::TCZ()]
                                  [DTI(CC::DIM())][DTI(CC::DIM())];

  for (unsigned int i = 0; i < CC::TCX(); ++i) {
    for (unsigned int j = 0; j < CC::TCY(); ++j) {
      for (unsigned int k = 0; k < CC::TCZ(); ++k) {
        for (unsigned int r = 0; r < DTI(CC::DIM()); ++r) {
          for (unsigned int c = 0; c < DTI(CC::DIM()); ++c) {
            velocity_gradient[i][j][k][r][c] = 0.0;
            shear_rate_tensor[i][j][k][r][c] = 0.0;
            shear_rate_tensor_squared[i][j][k][r][c] = 0.0;
          } // s
        }   // r
      }     // k
    }       // j
  }         // i

  // Calculate the gradient of the velocity vector
  BO::Stencils::ComputeVectorGradientAtCellCenter<DerivativeStencil>(
      u, v, w, cell_size, velocity_gradient);
  // compute tensor and squared tensor of the shear_rate
  ComputeShearRateTensor(velocity_gradient, shear_rate_tensor);
  TensorOperations::ComputeSquaredOfTensor(shear_rate_tensor,
                                           shear_rate_tensor_squared);

  // Computation of the actual shear rate
  for (unsigned int i = 0; i < CC::TCX(); ++i) {
    for (unsigned int j = 0; j < CC::TCY(); ++j) {
      for (unsigned int k = 0; k < CC::TCZ(); ++k) {

        double const second_invariant_of_shear_rate_tensor =
            TensorOperations::ComputeSecondInvariant(
                shear_rate_tensor, shear_rate_tensor_squared, i, j, k);

        // Formula: gamma_dot = 2*sqrt(I_2(S_ij))
        shear_rate[i][j][k] =
            2.0 * std::sqrt(std::abs(second_invariant_of_shear_rate_tensor));
      }
    }
  }
}

/**
 * @brief Computes the shear rate for a given velocity vectorial velocity field
 * at a certain position.
 * @param u,v,w Velocity buffer for which th shear rate is computed.
 * @param cell_size size of the cell for gradient computation.
 * @param i, j ,k indices at which position the shear rate should be computed.
 * @return shear rate at the given position.
 */
template <typename DerivativeStencil>
inline double
ComputeShearRate(double const (&u)[CC::TCX()][CC::TCY()][CC::TCZ()],
                 double const (&v)[CC::TCX()][CC::TCY()][CC::TCZ()],
                 double const (&w)[CC::TCX()][CC::TCY()][CC::TCZ()],
                 double const cell_size, unsigned int const i,
                 unsigned int const j, unsigned int const k) {

  // Calculate the velocity gradient
  std::array<std::array<double, 3>, 3> const velocity_gradient =
      SU::JacobianMatrix<DerivativeStencil>(u, v, w, i, j, k, cell_size);
  std::array<std::array<double, 3>, 3> shear_rate_tensor = {
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
  std::array<std::array<double, 3>, 3> shear_rate_tensor_squared = {
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};

  // compute tensor and squared tensor of the shear_rate
  ComputeShearRateTensor(velocity_gradient, shear_rate_tensor);
  TensorOperations::ComputeSquaredOfTensor(shear_rate_tensor,
                                           shear_rate_tensor_squared);

  // compute the second invariant of the shear rate tensor
  double const second_invariant_of_shear_rate_tensor =
      TensorOperations::ComputeSecondInvariant(shear_rate_tensor,
                                               shear_rate_tensor_squared);

  // retrun the shear_rate
  // Formula: gamma_dot = 2*sqrt(I_2(S_ij))
  return 2.0 * std::sqrt(std::abs(second_invariant_of_shear_rate_tensor));
}

} // namespace ShearRateOperations

#endif // SHEAR_RATE_UTILITIES_H
