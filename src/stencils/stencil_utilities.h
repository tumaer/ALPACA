//===----------------------- stencil_utilities.h --------------------------===//
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
#ifndef STENCIL_UTILITIES_H
#define STENCIL_UTILITIES_H

#include "apply_wrapper.h"

// Static assertions for the application of the reconstruction and derivative
// stencils
static_assert(ReconstructionStencilSetup::Concretize<
                  reconstruction_stencil>::type::DownstreamStencilSize() <
                  CC::HS(),
              "Halo size not enough for RECONSTRUCTION_STENCIL. Increase the "
              "halo size in compile_time_constants.h!");
static_assert(DerivativeStencilSetup::Concretize<
                  derivative_stencil>::type::DownstreamStencilSize() < CC::HS(),
              "Halo size not enough for DERIVATIVE_STENCIL. Increase the halo "
              "size in compile_time_constants.h!");

namespace StencilUtilities {

/**
 * @brief Applies the reconstruction operation of a given stencil for a linear
 * array.
 * @param array The array holding the information required for the stencil
 * evaluation.
 * @param cell_size The cell size used for stencil application.
 * @return The reconstructed value.
 * @tparam S The used stencil.
 * @tparam P The manner in which the stencil is applied.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, StencilProperty P, typename T = double,
          Dimension DIM = CC::DIM()>
constexpr T Reconstruction(std::array<double, S::StencilSize()> const &array,
                           double const cell_size) {
  return ApplyUtilities::Apply<S, P, T>(array, cell_size);
}

/**
 * @brief Applies the reconstruction operation of a given stencil for a single
 * cell of the given buffer.
 * @param buffer The buffer holding the information of all cells.
 * @param i,j,k The indices in the Cartesian directions for the cell that should
 * be evaluated.
 * @param cell_size The cell size used for stencil application.
 * @return The reconstructed value.
 * @tparam S The used stencil.
 * @tparam P The manner in which the stencil is applied.
 * @tparam D The direction in which the buffer should be applied.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, StencilProperty P, Direction D, typename T = double,
          Dimension DIM = CC::DIM()>
constexpr T Reconstruction(T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                           unsigned int const i, unsigned int const j,
                           unsigned int const k, double const cell_size) {
  return ApplyUtilities::Apply<S, P, D>(buffer, i, j, k, cell_size);
}

/**
 * @brief Applies the reconstruction operation of a given stencil for a linear
 * array in upwind direction.
 * @param array The array holding the information required for the stencil
 * evaluation.
 * @param upwind_decision Positive for UpwindLeft evaluation and negative for
 * UpwindRight.
 * @param cell_size The cell size used for stencil application.
 * @return The reconstructed value.
 * @tparam S The used stencil.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, typename T = double, Dimension DIM = CC::DIM()>
constexpr T
ReconstructionWithUpwinding(std::array<double, S::StencilSize()> const &array,
                            double const upwind_decision,
                            double const cell_size) {
  if (upwind_decision >= 0) {
    return ApplyUtilities::Apply<S, SP::UpwindLeft>(array, cell_size);
  } else {
    return ApplyUtilities::Apply<S, SP::UpwindRight>(array, cell_size);
  }
}

/**
 * @brief Applies the reconstruction operation of a given stencil for a single
 * cell of the given buffer in upwind direction.
 * @param buffer The buffer holding the information of all cells.
 * @param i,j,k The indices in the Cartesian directions for the cell that should
 * be evaluated.
 * @param upwind_decision Positive for UpwindLeft evaluation and negative for
 * UpwindRight.
 * @param cell_size The cell size used for stencil application.
 * @return The reconstructed value.
 * @tparam S The used stencil.
 * @tparam D The direction in which the buffer should be applied.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, Direction D, typename T = double,
          Dimension DIM = CC::DIM()>
constexpr T
ReconstructionWithUpwinding(T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                            unsigned int const i, unsigned int const j,
                            unsigned int const k, double const upwind_decision,
                            double const cell_size) {
  if (upwind_decision >= 0) {
    return ApplyUtilities::Apply<S, SP::UpwindLeft, D>(buffer, i, j, k,
                                                       cell_size);
  } else {
    return ApplyUtilities::Apply<S, SP::UpwindRight, D>(buffer, i, j, k,
                                                        cell_size);
  }
}

/**
 * @brief Computes the derivative of a given stencil for a linear array. The
 * stencils are evaluated central.
 * @param array The array holding the information required for the stencil
 * evaluation.
 * @param cell_size The cell size used for stencil application.
 * @return The value for the derivative.
 * @tparam S The used stencil.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, typename T = double, Dimension DIM = CC::DIM()>
constexpr T Derivative(std::array<double, S::StencilSize()> const &array,
                       double const cell_size) {
  return ApplyUtilities::Apply<S, SP::Central, T>(array, cell_size);
}

/**
 * @brief Computes the derivative of a given stencil for a single cell of the
 * given buffer. Derivative stencils are evaluated central. Reconstruction
 *        stencils as the mean of UpwindLeft and UpwindRight with a HjWeno
 * Scheme.
 * @param buffer The buffer holding the information of all cells.
 * @param i,j,k The indices in the Cartesian directions for the cell that should
 * be evaluated.
 * @param cell_size The cell size used for stencil application.
 * @return The reconstructed value.
 * @tparam S The used stencil.
 * @tparam D The direction in which the buffer should be applied.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, Direction D, typename T = double,
          Dimension DIM = CC::DIM()>
constexpr T Derivative(T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                       unsigned int const i, unsigned int const j,
                       unsigned int const k, double const cell_size) {
  switch (S::GetStencilType()) {
  case StencilType::Derivative:
    return ApplyUtilities::Apply<S, SP::Central, D>(buffer, i, j, k, cell_size);
  default:
    return 0.5 * (ApplyUtilities::ApplyHjWeno<S, SP::UpwindLeft, D>(
                      buffer, i, j, k, cell_size) +
                  ApplyUtilities::ApplyHjWeno<S, SP::UpwindRight, D>(
                      buffer, i, j, k, cell_size));
  }
}

/**
 * @brief Computes the derivative of a given stencil for a linear array. The
 * stencil property can be specified individually.
 * @param array The array holding the information required for the stencil
 * evaluation.
 * @param cell_size The cell size used for stencil application.
 * @return The value for the derivative.
 * @tparam S The used stencil.
 * @tparam P The manner in which the stencil is applied.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, StencilProperty P, typename T = double,
          Dimension DIM = CC::DIM()>
constexpr T Derivative(std::array<double, S::StencilSize()> const &array,
                       double const cell_size) {
  return ApplyUtilities::Apply<S, P, T>(array, cell_size);
}

/**
 * @brief Computes the derivative of a given stencil for a single cell of the
 * given buffer. It is differed between Derivative and Reconstruction stencils.
 *        The stencil property can be specified individually.
 * @param buffer The buffer holding the information of all cells.
 * @param i,j,k The indices in the Cartesian directions for the cell that should
 * be evaluated.
 * @param cell_size The cell size used for stencil application.
 * @return The reconstructed value.
 * @tparam S The used stencil.
 * @tparam P The manner in which the stencil is applied.
 * @tparam D The direction in which the buffer should be applied.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, StencilProperty P, Direction D, typename T = double,
          Dimension DIM = CC::DIM()>
constexpr T Derivative(T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                       unsigned int const i, unsigned int const j,
                       unsigned int const k, double const cell_size) {
  switch (S::GetStencilType()) {
  case StencilType::Derivative:
    return ApplyUtilities::Apply<S, P, D>(buffer, i, j, k, cell_size);
  default:
    return ApplyUtilities::ApplyHjWeno<S, P, D>(buffer, i, j, k, cell_size);
  }
}

/**
 * @brief Computes the derivative of a given stencil for a linear array in
 * upwind direction.
 * @param array The array holding the information required for the stencil
 * evaluation.
 * @param upwind_decision Positive for UpwindLeft evaluation and negative for
 * UpwindRight.
 * @param cell_size The cell size used for stencil application.
 * @return The value for the derivative.
 * @tparam S The used stencil.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, typename T = double, Dimension DIM = CC::DIM()>
constexpr T
DerivativeWithUpwinding(std::array<double, S::StencilSize()> const &array,
                        double const upwind_decision, double const cell_size) {
  if (upwind_decision >= 0) {
    return ApplyUtilities::Apply<S, SP::UpwindLeft, T>(array, cell_size);
  } else {
    return ApplyUtilities::Apply<S, SP::UpwindRight, T>(array, cell_size);
  }
}

/**
 * @brief Computes the derivative of a given stencil for a single cell of the
 * given buffer in upwind direction. It is differed between Derivative and
 * Reconstruction stencils.
 * @param buffer The buffer holding the information of all cells.
 * @param i,j,k The indices in the Cartesian directions for the cell that should
 * be evaluated.
 * @param cell_size The cell size used for stencil application.
 * @return The reconstructed value.
 * @tparam S The used stencil.
 * @tparam D The direction in which the buffer should be applied.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, Direction D, typename T = double,
          Dimension DIM = CC::DIM()>
constexpr T
DerivativeWithUpwinding(T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                        unsigned int const i, unsigned int const j,
                        unsigned int const k, double const upwind_decision,
                        double const cell_size) {
  switch (S::GetStencilType()) {
  case StencilType::Derivative: {
    if (upwind_decision >= 0) {
      return ApplyUtilities::Apply<S, SP::UpwindLeft, D>(buffer, i, j, k,
                                                         cell_size);
    } else {
      return ApplyUtilities::Apply<S, SP::UpwindRight, D>(buffer, i, j, k,
                                                          cell_size);
    }
  }
  default: {
    if (upwind_decision >= 0) {
      return ApplyUtilities::ApplyHjWeno<S, SP::UpwindLeft, D>(buffer, i, j, k,
                                                               cell_size);
    } else {
      return ApplyUtilities::ApplyHjWeno<S, SP::UpwindRight, D>(buffer, i, j, k,
                                                                cell_size);
    }
  }
  }
}

/**
 * @brief Computes the gradient (derivative in all Cartesian Directions) of a
 * given stencil for a single cell of the given buffer.
 * @param buffer The buffer holding the information of all cells.
 * @param i,j,k The indices in the Cartesian directions for the cell that should
 * be evaluated.
 * @param cell_size The cell size used for stencil application.
 * @return The gradient vector (3x1).
 * @tparam S The used stencil.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, typename T = double, Dimension DIM = CC::DIM()>
constexpr std::array<T, 3>
GradientVector(T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
               unsigned int const i, unsigned int const j, unsigned int const k,
               double const cell_size) {
  return {Derivative<S, Direction::X, T, DIM>(buffer, i, j, k, cell_size),
          DIM != Dimension::One
              ? Derivative<S, Direction::Y, T, DIM>(buffer, i, j, k, cell_size)
              : 0.0,
          DIM == Dimension::Three
              ? Derivative<S, Direction::Z, T, DIM>(buffer, i, j, k, cell_size)
              : 0.0};
}

/**
 * @brief Computes the jacobian matrix (derivative in all Cartesian Directions)
 * of a vectorial quantity of a given stencil for a single cell of the given
 * buffer.
 * @param buffer The buffer holding the information of all cells.
 * @param i,j,k The indices in the Cartesian directions for the cell that should
 * be evaluated.
 * @param cell_size The cell size used for stencil application.
 * @return The jacobian matrix (3x3).
 * @tparam S The used stencil.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, typename T = double, Dimension DIM = CC::DIM()>
constexpr std::array<std::array<T, 3>, 3>
JacobianMatrix(T const (&buffer_x)[CC::TCX()][CC::TCY()][CC::TCZ()],
               T const (&buffer_y)[CC::TCX()][CC::TCY()][CC::TCZ()],
               T const (&buffer_z)[CC::TCX()][CC::TCY()][CC::TCZ()],
               unsigned int const i, unsigned int const j, unsigned int const k,
               double const cell_size) {

  return {GradientVector<S, T, DIM>(buffer_x, i, j, k, cell_size),
          DIM != Dimension::One
              ? GradientVector<S, T, DIM>(buffer_y, i, j, k, cell_size)
              : std::array<T, 3>({0.0, 0.0, 0.0}),
          DIM == Dimension::Three
              ? GradientVector<S, T, DIM>(buffer_z, i, j, k, cell_size)
              : std::array<T, 3>({0.0, 0.0, 0.0})};
}

/**
 * @brief Computes the curl of a vectorial quantity of a given stencil for a
 * single cell of the given buffer.
 * @param buffer The buffer holding the information of all cells.
 * @param i,j,k The indices in the Cartesian directions for the cell that should
 * be evaluated.
 * @param cell_size The cell size used for stencil application.
 * @return The curl (3x1).
 * @tparam S The used stencil.
 * @tparam T The type of the output value (default: double).
 * @tparam DIM The dimension the stencil is applied on (default: CC::DIM())
 */
template <typename S, typename T = double, Dimension DIM = CC::DIM()>
constexpr std::array<T, 3>
Curl(T const (&buffer_x)[CC::TCX()][CC::TCY()][CC::TCZ()],
     T const (&buffer_y)[CC::TCX()][CC::TCY()][CC::TCZ()],
     T const (&buffer_z)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
     unsigned int const j, unsigned int const k, double const cell_size) {
  std::array<std::array<T, 3>, 3> const gradient = JacobianMatrix<S, T, DIM>(
      buffer_x, buffer_y, buffer_z, i, j, k, cell_size);

  return {gradient[2][1] - gradient[1][2], gradient[0][2] - gradient[2][0],
          gradient[1][0] - gradient[0][1]};
}

} // namespace StencilUtilities

namespace SU = StencilUtilities;

#endif // STENCIL_UTILITIES_H
