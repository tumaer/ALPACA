//===----------------- state_reconstruction_settings.h --------------------===//
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
#ifndef STATE_RECONSTRUCTION_SETTINGS_H
#define STATE_RECONSTRUCTION_SETTINGS_H

enum class StateReconstructionType { Conservative, Primitive, Characteristic };

constexpr StateReconstructionType state_reconstruction_type =
    StateReconstructionType::Characteristic;

/**
 * @brief Provides a string representation of the state reconstruction type.
 * @param reconstruction State reconstruction set to be stringified.
 * @return String of the given state reconstruction type.
 */
inline std::string SetToString(StateReconstructionType const reconstruction) {
  switch (reconstruction) {
  case StateReconstructionType::Conservative:
    return "Conservative";
  case StateReconstructionType::Primitive:
    return "Primitive";
  case StateReconstructionType::Characteristic:
    return "Characteristic";
  default:
    return "ERROR: This reconstruction type is not (yet) defined!";
  }
}

#endif // STATE_RECONSTRUCTION_SETTINGS_H
