//===----------------------- geometry_settings.h --------------------------===//
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
#ifndef GEOMETRY_SETTINGS_H
#define GEOMETRY_SETTINGS_H

/**
 * @brief Identifer for the cut cell criteria.
 */
enum class CutCellCriteria { SignChangeBased, ValueBased };

/**
 * @brief Identifer whther to use differentiation or reconstruction stencil
 * inside the geometry calculator.
 */
enum class GeometryStencilType { Reconstruction, Derivative };

#endif // GEOMETRY_SETTINGS_H
