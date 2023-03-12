//===------------------- state_reconstruction_setup.h ---------------------===//
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
#ifndef STATE_RECONSTRUCTION_SETUP_H
#define STATE_RECONSTRUCTION_SETUP_H

#include "solvers/state_reconstruction/characteristic_state_reconstruction.h"
#include "solvers/state_reconstruction/conservative_state_reconstruction.h"
#include "solvers/state_reconstruction/gamma_characteristic_state_reconstruction.h"
#include "solvers/state_reconstruction/gamma_primitive_state_reconstruction.h"
#include "solvers/state_reconstruction/primitive_state_reconstruction.h"
#include "user_specifications/equation_settings.h"
#include "user_specifications/state_reconstruction_settings.h"

static_assert(
    !(active_equations == EquationSet::Isentropic &&
      state_reconstruction_type == StateReconstructionType::Characteristic),
    "Characteristic reconstruction not implemented for isentropic equations!");

/**
 * @brief A namespace to get a StateReconstruction type based on a specified
 * constexpr.
 */
namespace StateReconstructionSetup {

namespace EulerNavierStokes {
/**
 * @brief Function returning the Euler or Navier-Stokes equations state
 * reconstruction matching the type in the template argument.
 * @tparam StateReconstructionType Specification of the state reconstruction
 * type.
 */
template <StateReconstructionType> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<StateReconstructionType::Conservative> {
  using type = ConservativeStateReconstruction;
};

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<StateReconstructionType::Primitive> {
  using type = PrimitiveStateReconstruction;
};
/**
 * @brief See generic implementation.
 */
template <> struct Concretize<StateReconstructionType::Characteristic> {
  using type = CharacteristicStateReconstruction;
};
} // namespace EulerNavierStokes

namespace Isentropic {
/**
 * @brief Function returning the Isentropic state reconstruction matching the
 * type in the template argument.
 * @tparam StateReconstructionType Specification of the state reconstruction
 * type.
 */
template <StateReconstructionType> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<StateReconstructionType::Conservative> {
  using type = ConservativeStateReconstruction;
};

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<StateReconstructionType::Primitive> {
  using type = PrimitiveStateReconstruction;
};
} // namespace Isentropic

namespace GammaModel {
/**
 * @brief Function returning the Gamma Model state reconstruction matching the
 * type in the template argument.
 * @tparam StateReconstructionType Specification of the state reconstruction
 * type.
 */
template <StateReconstructionType> struct Concretize;

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<StateReconstructionType::Conservative> {
  using type = ConservativeStateReconstruction;
};

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<StateReconstructionType::Primitive> {
  using type = GammaPrimitiveStateReconstruction;
};

/**
 * @brief See generic implementation.
 */
template <> struct Concretize<StateReconstructionType::Characteristic> {
  using type = GammaCharacteristicStateReconstruction;
};
} // namespace GammaModel

/**
 * @brief Function returning the state reconstruction type matching the type in
 * the template argument according to the equation set in the second template
 * parameter.
 * @tparam StateReconstructionType Specification of the state reconstruction
 * type.
 * @tparam EquationSet Specification of the equation(s) beeing solved.
 */
template <StateReconstructionType S, EquationSet> struct Dispatch {
  using type = typename EulerNavierStokes::Concretize<S>::type;
};

/**
 * @brief See generic implementation.
 */
template <StateReconstructionType S>
struct Dispatch<S, EquationSet::Isentropic> {
  using type = typename Isentropic::Concretize<S>::type;
};

/**
 * @brief See generic implementation.
 */
template <StateReconstructionType S>
struct Dispatch<S, EquationSet::GammaModel> {
  using type = typename GammaModel::Concretize<S>::type;
};

/**
 * @brief Function returning the state reconstruction type for the (globally)
 * selected equation and the given reconstruction template argument.
 * @tparam StateReconstructionType Specification of the state reconstruction
 * type.
 */
template <StateReconstructionType S> struct Concretize {
  using type = typename Dispatch<S, active_equations>::type;
};

} // namespace StateReconstructionSetup

#endif // STATE_RECONSTRUCTION_SETUP_H
