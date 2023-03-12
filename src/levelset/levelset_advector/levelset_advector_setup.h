//===------------------- levelset_advector_setup.h ------------------------===//
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
#ifndef LEVELSET_ADVECTOR_SETUP_H
#define LEVELSET_ADVECTOR_SETUP_H

#include "derivative_stencil_single_levelset_advector.h"
#include "hj_derivative_stencil_single_levelset_advector.h"
#include "hj_reconstruction_stencil_single_levelset_advector.h"
#include "reconstruction_stencil_single_levelset_advector.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief A namespace to get a LevelsetAdvector type based on a specified
 * constexpr.
 */
namespace LevelsetAdvectorSetup {

/**
 * @brief Function returning the typedef of a LevelsetAdvector based on a
 * constexpr template.
 *
 * @tparam LevelsetAdvectors The constexpr template parameter to specify the
 * exact LevelsetAdvector type.
 */
template <LevelsetAdvectors> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<LevelsetAdvectors::DerivativeStencil> {
  typedef DerivativeStencilSingleLevelsetAdvector type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<LevelsetAdvectors::ReconstructionStencil> {
  typedef ReconstructionStencilSingleLevelsetAdvector type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<LevelsetAdvectors::HjReconstructionStencil> {
  typedef HjReconstructionStencilSingleLevelsetAdvector type;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<LevelsetAdvectors::HjDerivativeStencil> {
  typedef HjDerivativeStencilSingleLevelsetAdvector type;
};

} // namespace LevelsetAdvectorSetup

#endif // LEVELSET_ADVECTOR_SETUP_H
