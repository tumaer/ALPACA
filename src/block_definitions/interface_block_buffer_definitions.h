//===------------- interface_block_buffer_definitions.h -------------------===//
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
#ifndef INTERFACE_BLOCK_BUFFER_DEFINITIONS_H
#define INTERFACE_BLOCK_BUFFER_DEFINITIONS_H

#include "block_definitions/field_buffer.h"

/**
 * @brief Identifier of all buffers in an interface block (to allow single
 * buffer access, e.g., in halo update or extension).
 */
enum class InterfaceBlockBufferType {
  // interface description buffers
  LevelsetBase,
  LevelsetRightHandSide,
  LevelsetReinitialized,
  LevelsetIntegrated,
  VolumeFractionBase,
  VolumeFractionRightHandSide,
  VolumeFractionReinitialized,
  VolumeFractionIntegrated,
  // interface state buffers
  InterfaceStateVelocity,
  InterfaceStatePressurePositive,
  InterfaceStatePressureNegative,
  // interface parameter buffers
  InterfaceParameterSurfaceTensionCoefficient
};

/**
 * @brief Mapping function to map an interface description buffer to the given
 * interface block buffer type.
 * @param description_field The Identifier of the description field (Levelset or
 * VolumeFraction).
 * @param buffer_type The corresponding interface description buffer type (Base,
 * RightHandSide, Reinitialized).
 * @return Appropriate InterfaceBlockBufferType.
 */
static constexpr InterfaceBlockBufferType
MapInterfaceDescritpionToInterfaceBlockBufferType(
    InterfaceDescription const description_field,
    InterfaceDescriptionBufferType const buffer_type) {
  switch (buffer_type) {
  case InterfaceDescriptionBufferType::Base: {
    switch (description_field) {
    case InterfaceDescription::Levelset: {
      return InterfaceBlockBufferType::LevelsetBase;
    }
    default: {
      return InterfaceBlockBufferType::VolumeFractionBase;
    }
    }
  } break;
  case InterfaceDescriptionBufferType::RightHandSide: {
    switch (description_field) {
    case InterfaceDescription::Levelset: {
      return InterfaceBlockBufferType::LevelsetRightHandSide;
    }
    default: {
      return InterfaceBlockBufferType::VolumeFractionRightHandSide;
    }
    }
  } break;
  case InterfaceDescriptionBufferType::Integrated: {
    switch (description_field) {
    case InterfaceDescription::Levelset: {
      return InterfaceBlockBufferType::LevelsetIntegrated;
    }
    default: {
      return InterfaceBlockBufferType::VolumeFractionIntegrated;
    }
    }
  } break;
  // Last possibility (reinitialized)
  default: {
    switch (description_field) {
    case InterfaceDescription::Levelset: {
      return InterfaceBlockBufferType::LevelsetReinitialized;
    }
    default: {
      return InterfaceBlockBufferType::VolumeFractionReinitialized;
    }
    }
  }
  }
}

/**
 * @brief Mapping function to map an interface state to the given interface
 * block buffer type.
 * @param state The Identifier of the interface state field (Velocity,
 * PressurePositive or PressureNegative).
 * @return Appropriate InterfaceBlockBufferType.
 */
static constexpr InterfaceBlockBufferType
MapInterfaceStateToInterfaceBlockBufferType(InterfaceState const state) {
  switch (state) {
  case InterfaceState::Velocity: {
    return InterfaceBlockBufferType::InterfaceStateVelocity;
  }
  case InterfaceState::PressurePositive: {
    return InterfaceBlockBufferType::InterfaceStatePressurePositive;
  }
  // Last possibility
  default: {
    return InterfaceBlockBufferType::InterfaceStatePressureNegative;
  }
  }
}

/**
 * @brief Mapping function to map an interface parameter to the given interface
 * block buffer type.
 * @param parameter The Identifier of the interface parameter field
 * (SurfaceTensionCoefficient).
 * @return Appropriate InterfaceBlockBufferType.
 */
static constexpr InterfaceBlockBufferType
MapInterfaceParameterToInterfaceBlockBufferType(
    InterfaceParameter const parameter) {
  switch (parameter) {
  default: {
    return InterfaceBlockBufferType::
        InterfaceParameterSurfaceTensionCoefficient;
  }
  }
}

/**
 * @brief Mapping function to map an interface field and index to the given
 * interface block buffer type.
 * @param field_type The Identifier of the interface field (States, Parameters
 * or Description).
 * @param field_index Index of the given field.
 * @param buffer_type The corresponding interface buffer type (Base,
 * RightHandSide, Reinitialized).
 * @return Appropriate InterfaceBlockBufferType.
 */
static constexpr InterfaceBlockBufferType
MapInterfaceFieldToInterfaceBlockBufferType(
    InterfaceFieldType const field_type, unsigned int const field_index,
    InterfaceDescriptionBufferType const buffer_type =
        InterfaceDescriptionBufferType::RightHandSide) {
  switch (field_type) {
  case InterfaceFieldType::Description: {
    return MapInterfaceDescritpionToInterfaceBlockBufferType(
        IF::ITID(field_index), buffer_type);
  }
  case InterfaceFieldType::Parameters: {
    return MapInterfaceParameterToInterfaceBlockBufferType(
        IF::ITIP(field_index));
  }
  default: {
    return MapInterfaceStateToInterfaceBlockBufferType(IF::ITIS(field_index));
  }
  }
}

#endif // INTERFACE_BLOCK_BUFFER_DEFINITIONS_H
