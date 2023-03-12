//===------------------- internal_boundary_types.h ------------------------===//
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
/**
 * @brief Unique identifier to separate the kind of Internal boundary with
 * respect to its MPI treatment and topology information.
 */
enum class InternalBoundaryType {
  JumpBoundaryLocal,
  JumpBoundaryMpiSend,
  JumpBoundaryMpiRecv,
  NoJumpBoundaryLocal,
  NoJumpBoundaryMpiSend,
  NoJumpBoundaryMpiRecv
};
