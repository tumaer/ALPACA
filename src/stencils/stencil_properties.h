//===---------------------- stencil_properties.h --------------------------===//
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
#ifndef STENCIL_PROPERTIES_H
#define STENCIL_PROPERTIES_H

/**
 * @brief Unique identifier whether a stencil is a derivative or reconstruction
 * stencil.
 */
enum class StencilType { Reconstruction, Derivative };

/**
 * @brief Unique identifier to indicate whether a stencil should be applied
 * UpwindLeft, UpwindRight or Central.
 */
enum class StencilProperty { UpwindLeft, UpwindRight, Central };
using SP = StencilProperty;

/**
 * @brief Unique identifier to indicate whether a stencil should reconstruct
 * based on cell-centered values or based on differences (HJ-WENO like
 * reconstruction).
 */
enum class ReconstructionType { Values, Differences };

#endif // STENCIL_PROPERTIES_H
