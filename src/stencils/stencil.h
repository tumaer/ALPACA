//===--------------------------- stencil.h --------------------------------===//
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
#ifndef SPATIAL_RECONSTRUCTION_STENCIL_H
#define SPATIAL_RECONSTRUCTION_STENCIL_H

#include <limits>
#include <vector>

#include "enums/direction_definition.h"
#include "stencils/stencil_properties.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief The SpatialReconstructionStencil class provides an interface for
 * computational Stencils needed throughout the simulation.
 */
template <typename DerivedStencil> class Stencil {

  friend DerivedStencil;

  explicit constexpr Stencil() = default;

protected:
  static constexpr double epsilon_ = std::numeric_limits<double>::epsilon();

public:
  ~Stencil() = default;
  Stencil(Stencil const &) = delete;
  Stencil &operator=(Stencil const &) = delete;
  Stencil(Stencil &&) = delete;
  Stencil &operator=(Stencil &&) = delete;

  /**
   * @brief Gives the number of cells needed for a single stencil evaluation.
   * @return Size of the stencil, i.e. number of data cells the stencil works
   * on.
   */
  static constexpr unsigned int StencilSize() {
    return DerivedStencil::stencil_size_;
  }
  /**
   * @brief Gives the size of the stencil in down stream direction.
   * @return Size of the stencil arm reaching down stream, i.e. number of data
   * cells that lay down stream the stencil works on.
   */
  static constexpr unsigned int DownstreamStencilSize() {
    return DerivedStencil::downstream_stencil_size_;
  }
  /**
   * @brief Return the type of a stencil.
   * @return The type of a stencil.
   */
  static constexpr StencilType GetStencilType() {
    return DerivedStencil::stencil_type_;
  }

  /**
   * @brief Applies the SpatialReconstructionStencil to the provided Array.
   * @param array The array on which to apply the spatial reconstruction
   * stencil.
   * @param cell_size The cell size of the corresponding block.
   * @return Value at the position of interest.
   * @tparam S The used stencil.
   * @note Hotpath function.
   * @note The stencil template is used to instantiate the stencil on the fly.
   */
  template <typename S>
  constexpr double Apply(std::array<double, S::StencilSize()> const &array,
                         std::array<int const, 2> const evaluation_properties,
                         double const cell_size) const {
    return static_cast<DerivedStencil const &>(*this).ApplyImplementation(
        array, evaluation_properties, cell_size);
  }
};

#endif // SPATIAL_RECONSTRUCTION_STENCIL_H
