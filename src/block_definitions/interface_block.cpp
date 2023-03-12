//===----------------------- interface_block.cpp --------------------------===//
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
#include "block_definitions/interface_block.h"
#include "utilities/buffer_operations.h"
#include <stdexcept>

/**
 * @brief Constructor to create a interface block according to an already
 * computed levelset field.
 * @param levelset_initial Reference to array holding the levelset field to be
 * in this interface block.
 */
InterfaceBlock::InterfaceBlock(
    double const (&levelset_initial)[CC::TCX()][CC::TCY()][CC::TCZ()]) {

  // In the base buffer all values are set to zero
  BO::SetFieldBuffer(GetBaseBuffer(), 0.0);
  // In the right-hand and reinitialized buffer the levelset is copied and
  // volume fraction is set to zero
  BO::CopySingleBuffer(levelset_initial,
                       GetRightHandSideBuffer(InterfaceDescription::Levelset));
  BO::CopySingleBuffer(levelset_initial,
                       GetReinitializedBuffer(InterfaceDescription::Levelset));
  BO::CopySingleBuffer(levelset_initial,
                       GetIntegratedBuffer(InterfaceDescription::Levelset));
  BO::SetSingleBuffer(
      GetRightHandSideBuffer(InterfaceDescription::VolumeFraction), 0.0);
  BO::SetSingleBuffer(
      GetReinitializedBuffer(InterfaceDescription::VolumeFraction), 0.0);
  BO::SetSingleBuffer(GetIntegratedBuffer(InterfaceDescription::VolumeFraction),
                      0.0);
  // In the initial buffer all values are set to zero
  BO::SetFieldBuffer(GetInitialBuffer(), 0.0);
  // In the interface state buffer all values are set to zero
  BO::SetFieldBuffer(GetInterfaceStateBuffer(), 0.0);

  // Only set buffer of parameter to zero if they are present
  if constexpr (CC::InterfaceParameterModelActive()) {
    BO::SetFieldBuffer(GetInterfaceParameterBuffer(), 0.0);
  }
}

/**
 * @brief Constructor to create an initial homogenous levelset field on the
 * interface block.
 * @param levelset_initial The value to be imposed as levelset function.
 */
InterfaceBlock::InterfaceBlock(double const levelset_initial) {

  // In the base buffer only the volume fraction is set to uniform values, the
  // levelset is set to zero
  BO::SetSingleBuffer(GetBaseBuffer(InterfaceDescription::Levelset), 0.0);
  BO::SetSingleBuffer(GetBaseBuffer(InterfaceDescription::VolumeFraction),
                      levelset_initial > 0 ? 1.0 : 0.0);
  // In the right-hand and reinitialized buffer the levelset is set to uniform
  // given value and volume fraction is set to zero
  BO::SetSingleBuffer(GetRightHandSideBuffer(InterfaceDescription::Levelset),
                      levelset_initial);
  BO::SetSingleBuffer(GetReinitializedBuffer(InterfaceDescription::Levelset),
                      levelset_initial);
  BO::SetSingleBuffer(GetIntegratedBuffer(InterfaceDescription::Levelset),
                      levelset_initial);
  BO::SetSingleBuffer(
      GetRightHandSideBuffer(InterfaceDescription::VolumeFraction), 0.0);
  BO::SetSingleBuffer(
      GetReinitializedBuffer(InterfaceDescription::VolumeFraction),
      levelset_initial > 0 ? 1.0 : 0.0);
  BO::SetSingleBuffer(GetIntegratedBuffer(InterfaceDescription::VolumeFraction),
                      0.0);
  // In the initial buffer all values are set to zero
  BO::SetFieldBuffer(GetInitialBuffer(), 0.0);
  // In the interface state buffer all values are set to zero
  BO::SetFieldBuffer(GetInterfaceStateBuffer(), 0.0);

  // Only set buffer of parameter to zero if they are present
  if constexpr (CC::InterfaceParameterModelActive()) {
    BO::SetFieldBuffer(GetInterfaceParameterBuffer(), 0.0);
  }
}

/**
 * @brief Gives a reference to the corresponding buffer.
 * @param field_type The interface field type of the buffer.
 * @param field_index The index of the field asked for.
 * @param buffer_type If a interface field is wanted, the interface buffer type.
 * Defaults to InterfaceDescriptionBufferType::RightHandSide.
 * @return Reference to Array that is the requested buffer.
 */
auto InterfaceBlock::GetFieldBuffer(
    InterfaceFieldType const field_type, unsigned int const field_index,
    InterfaceDescriptionBufferType const buffer_type)
    -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  switch (field_type) {
  case InterfaceFieldType::Description: {
    return GetInterfaceDescriptionBuffer(buffer_type)[field_index];
  }
  case InterfaceFieldType::Parameters: {
    return parameters_[field_index];
  }
  default: { // InterfaceFieldType::States:
    return states_[field_index];
  }
  }
}

/**
 * @brief Const overload.
 */
auto InterfaceBlock::GetFieldBuffer(
    InterfaceFieldType const field_type, unsigned int const field_index,
    InterfaceDescriptionBufferType const buffer_type) const
    -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  switch (field_type) {
  case InterfaceFieldType::Description: {
    return GetInterfaceDescriptionBuffer(buffer_type)[field_index];
  }
  case InterfaceFieldType::Parameters: {
    return parameters_[field_index];
  }
  default: { // InterfaceFieldType::States:
    return states_[field_index];
  }
  }
}

/**
 * @brief Gives a reference to the corresponding Base buffer.
 * @param interface_description Decider which buffer is to be returned.
 * @return Reference to the array, that is the requested buffer.
 */
auto InterfaceBlock::GetBaseBuffer(
    InterfaceDescription const interface_description)
    -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return base_[interface_description];
}

/**
 * @brief Const overload.
 */
auto InterfaceBlock::GetBaseBuffer(
    InterfaceDescription const interface_description) const
    -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return base_[interface_description];
}

/**
 * @brief Gives a Reference to the corresponding right-hand side buffer.
 * @param interface_description Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto InterfaceBlock::GetRightHandSideBuffer(
    InterfaceDescription const interface_description)
    -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return right_hand_side_[interface_description];
}

/**
 * @brief Const overload.
 */
auto InterfaceBlock::GetRightHandSideBuffer(
    InterfaceDescription const interface_description) const
    -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return right_hand_side_[interface_description];
}

/**
 * @brief Gives a Reference to the corresponding reinitialized buffer.
 * @param interface_description Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto InterfaceBlock::GetReinitializedBuffer(
    InterfaceDescription const interface_description)
    -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return reinitialized_[interface_description];
}

/**
 * @brief Const overload.
 */
auto InterfaceBlock::GetReinitializedBuffer(
    InterfaceDescription const interface_description) const
    -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return reinitialized_[interface_description];
}

/**
 * @brief Gives a Reference to the corresponding initial buffer.
 * @param interface_description Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto InterfaceBlock::GetInitialBuffer(
    InterfaceDescription const interface_description)
    -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return initial_[interface_description];
}

/**
 * @brief Const overload.
 */
auto InterfaceBlock::GetInitialBuffer(
    InterfaceDescription const interface_description) const
    -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return initial_[interface_description];
}

/**
 * @brief Gives a Reference to the corresponding integrated buffer.
 * @param interface_description Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto InterfaceBlock::GetIntegratedBuffer(
    InterfaceDescription const interface_description)
    -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return integrated_[interface_description];
}

/**
 * @brief Const overload.
 */
auto InterfaceBlock::GetIntegratedBuffer(
    InterfaceDescription const interface_description) const
    -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return integrated_[interface_description];
}

/**
 * @brief Wrapper function that returns an interface description buffer (base,
 * right-hand side, reinitialized, initial, or integrated). The decision is made
 * based on the template parameter. Implementation for the base buffer.
 * @return Base interface description struct.
 */
template <>
InterfaceDescriptions &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::Base>() {
  return GetBaseBuffer();
}

/**
 * @brief Const overload.
 */
template <>
InterfaceDescriptions const &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::Base>() const {
  return GetBaseBuffer();
}

/**
 * @brief Wrapper function that returns an interface description buffer (base,
 * right-hand side, reinitialized, initial, or integrated). The decision is made
 * based on the template parameter. Implementation for the right-hand side
 * buffer.
 * @return Right-hand side interface description struct.
 */
template <>
InterfaceDescriptions &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::RightHandSide>() {
  return GetRightHandSideBuffer();
}

/**
 * @brief Const overload.
 */
template <>
InterfaceDescriptions const &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::RightHandSide>() const {
  return GetRightHandSideBuffer();
}

/**
 * @brief Wrapper function that returns an interface description buffer (base,
 * right-hand side, reinitialized, initial, or integrated). The decision is made
 * based on the template parameter. Implementation for the reinitialized buffer.
 * @return Reinitialized interface description struct.
 */
template <>
InterfaceDescriptions &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::Reinitialized>() {
  return GetReinitializedBuffer();
}

/**
 * @brief Const overload.
 */
template <>
InterfaceDescriptions const &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::Reinitialized>() const {
  return GetReinitializedBuffer();
}

/**
 * @brief Wrapper function that returns an interface description buffer (base,
 * right-hand side, reinitialized, initial, or integrated). The decision is made
 * based on the template parameter. Implementation for the initial buffer.
 * @return Initial interface description struct.
 */
template <>
InterfaceDescriptions &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::Initial>() {
  return GetInitialBuffer();
}

/**
 * @brief Const overload.
 */
template <>
InterfaceDescriptions const &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::Initial>() const {
  return GetInitialBuffer();
}

/**
 * @brief Wrapper function that returns an interface description buffer (base,
 * right-hand side, reinitialized, initial, or integrated). The decision is made
 * based on the template parameter. Implementation for the integrated buffer.
 * @return Integrated interface description struct.
 */
template <>
InterfaceDescriptions &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::Integrated>() {
  return GetIntegratedBuffer();
}

/**
 * @brief Const overload.
 */
template <>
InterfaceDescriptions const &InterfaceBlock::GetInterfaceDescriptionBuffer<
    InterfaceDescriptionBufferType::Integrated>() const {
  return GetIntegratedBuffer();
}

/**
 * @brief Gives access to the base buffer.
 * @return base buffer struct.
 */
InterfaceDescriptions &InterfaceBlock::GetBaseBuffer() { return base_; }

/**
 * @brief Const overload.
 */
InterfaceDescriptions const &InterfaceBlock::GetBaseBuffer() const {
  return base_;
}

/**
 * @brief Gives access to the right-hand side buffer.
 * @return Right hand side buffer struct.
 */
InterfaceDescriptions &InterfaceBlock::GetRightHandSideBuffer() {
  return right_hand_side_;
}

/**
 * @brief Const overload.
 */
InterfaceDescriptions const &InterfaceBlock::GetRightHandSideBuffer() const {
  return right_hand_side_;
}

/**
 * @brief Gives access to the reinitialized buffer.
 * @return Reinitialized buffer struct.
 */
InterfaceDescriptions &InterfaceBlock::GetReinitializedBuffer() {
  return reinitialized_;
}

/**
 * @brief Const overload.
 */
InterfaceDescriptions const &InterfaceBlock::GetReinitializedBuffer() const {
  return reinitialized_;
}

/**
 * @brief Gives access to the initial buffer.
 * @return initial buffer struct.
 */
InterfaceDescriptions &InterfaceBlock::GetInitialBuffer() { return initial_; }

/**
 * @brief Const overload.
 */
InterfaceDescriptions const &InterfaceBlock::GetInitialBuffer() const {
  return initial_;
}

/**
 * @brief Gives access to the integrated buffer.
 * @return integrated buffer struct.
 */
InterfaceDescriptions &InterfaceBlock::GetIntegratedBuffer() {
  return integrated_;
}

/**
 * @brief Const overload.
 */
InterfaceDescriptions const &InterfaceBlock::GetIntegratedBuffer() const {
  return integrated_;
}

/**
 * @brief Gives access to the interface descritpion buffer of given type.
 * @param buffer_type InterfaceDescription type of the buffer asked for.
 * @return buffer struct of given type.
 */
InterfaceDescriptions &InterfaceBlock::GetInterfaceDescriptionBuffer(
    InterfaceDescriptionBufferType const buffer_type) {
  switch (buffer_type) {
  case InterfaceDescriptionBufferType::RightHandSide: {
    return right_hand_side_;
  }
  case InterfaceDescriptionBufferType::Reinitialized: {
    return reinitialized_;
  }
  case InterfaceDescriptionBufferType::Base: {
    return base_;
  }
  case InterfaceDescriptionBufferType::Integrated: {
    return integrated_;
  }
  default: { // case InterfaceDescriptionBufferType::Initial:
    return initial_;
  }
  }
}

/**
 * @brief Const overload.
 */
InterfaceDescriptions const &InterfaceBlock::GetInterfaceDescriptionBuffer(
    InterfaceDescriptionBufferType const buffer_type) const {
  switch (buffer_type) {
  case InterfaceDescriptionBufferType::RightHandSide: {
    return right_hand_side_;
  }
  case InterfaceDescriptionBufferType::Reinitialized: {
    return reinitialized_;
  }
  case InterfaceDescriptionBufferType::Base: {
    return base_;
  }
  case InterfaceDescriptionBufferType::Integrated: {
    return integrated_;
  }
  default: { // case InterfaceDescriptionBufferType::Initial:
    return initial_;
  }
  }
}

/**
 * @brief Gives a Reference to the corresponding interface state buffer.
 * @param state_type Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto InterfaceBlock::GetInterfaceStateBuffer(InterfaceState const state_type)
    -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return states_[state_type];
}

/**
 * @brief Const overload.
 */
auto InterfaceBlock::GetInterfaceStateBuffer(InterfaceState const state_type)
    const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return states_[state_type];
}

/**
 * @brief Gives access to the interface state struct.
 * @return Interface states buffer struct.
 */
InterfaceStates &InterfaceBlock::GetInterfaceStateBuffer() { return states_; }

/**
 * @brief Const overload.
 */
InterfaceStates const &InterfaceBlock::GetInterfaceStateBuffer() const {
  return states_;
}

/**
 * @brief Gives a Reference to the corresponding interface Parameter buffer.
 * @param parameter_type Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto InterfaceBlock::GetInterfaceParameterBuffer(
    InterfaceParameter const parameter_type)
    -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return parameters_[parameter_type];
}

/**
 * @brief Const overload.
 */
auto InterfaceBlock::GetInterfaceParameterBuffer(
    InterfaceParameter const parameter_type) const
    -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  return parameters_[parameter_type];
}

/**
 * @brief Gives access to the interface parameters struct.
 * @return Interface parameter buffer struct.
 */
InterfaceParameters &InterfaceBlock::GetInterfaceParameterBuffer() {
  return parameters_;
}

/**
 * @brief Const overload.
 */
InterfaceParameters const &InterfaceBlock::GetInterfaceParameterBuffer() const {
  return parameters_;
}

/**
 * @brief Gives the requested buffer of a specific single interface block
 * buffer.
 * @param buffer_type type of the interface block buffer requested.
 */
auto InterfaceBlock::GetBuffer(InterfaceBlockBufferType const buffer_type)
    -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  switch (buffer_type) {
  // interface descriptions
  case InterfaceBlockBufferType::LevelsetBase: {
    return base_[InterfaceDescription::Levelset];
  }
  case InterfaceBlockBufferType::VolumeFractionBase: {
    return base_[InterfaceDescription::VolumeFraction];
  }
  case InterfaceBlockBufferType::LevelsetRightHandSide: {
    return right_hand_side_[InterfaceDescription::Levelset];
  }
  case InterfaceBlockBufferType::VolumeFractionRightHandSide: {
    return right_hand_side_[InterfaceDescription::VolumeFraction];
  }
  case InterfaceBlockBufferType::LevelsetReinitialized: {
    return reinitialized_[InterfaceDescription::Levelset];
  }
  case InterfaceBlockBufferType::VolumeFractionReinitialized: {
    return reinitialized_[InterfaceDescription::VolumeFraction];
  }
  case InterfaceBlockBufferType::LevelsetIntegrated: {
    return integrated_[InterfaceDescription::Levelset];
  }
  case InterfaceBlockBufferType::VolumeFractionIntegrated: {
    return integrated_[InterfaceDescription::VolumeFraction];
  }
  // interface states
  case InterfaceBlockBufferType::InterfaceStateVelocity: {
    return states_[InterfaceState::Velocity];
  }
  case InterfaceBlockBufferType::InterfaceStatePressurePositive: {
    return states_[InterfaceState::PressurePositive];
  }
  case InterfaceBlockBufferType::InterfaceStatePressureNegative: {
    return states_[InterfaceState::PressureNegative];
  }
  // interface parameters
#ifdef PERFORMANCE
  default: {
    return parameters_[InterfaceParameter::SurfaceTensionCoefficient];
  }
#else
  case InterfaceBlockBufferType::InterfaceParameterSurfaceTensionCoefficient: {
    return parameters_[InterfaceParameter::SurfaceTensionCoefficient];
  }
  default: {
    throw std::invalid_argument(
        "Requested buffer does not exist (impossible error");
  } break;
#endif
  }
}

/**
 * @brief Const overload.
 */
auto InterfaceBlock::GetBuffer(InterfaceBlockBufferType const buffer_type) const
    -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  switch (buffer_type) {
  // interface description cases
  case InterfaceBlockBufferType::LevelsetBase: {
    return base_[InterfaceDescription::Levelset];
  }
  case InterfaceBlockBufferType::VolumeFractionBase: {
    return base_[InterfaceDescription::VolumeFraction];
  }
  case InterfaceBlockBufferType::LevelsetRightHandSide: {
    return right_hand_side_[InterfaceDescription::Levelset];
  }
  case InterfaceBlockBufferType::VolumeFractionRightHandSide: {
    return right_hand_side_[InterfaceDescription::VolumeFraction];
  }
  case InterfaceBlockBufferType::LevelsetReinitialized: {
    return reinitialized_[InterfaceDescription::Levelset];
  }
  case InterfaceBlockBufferType::VolumeFractionReinitialized: {
    return reinitialized_[InterfaceDescription::VolumeFraction];
  }
  case InterfaceBlockBufferType::LevelsetIntegrated: {
    return integrated_[InterfaceDescription::Levelset];
  }
  case InterfaceBlockBufferType::VolumeFractionIntegrated: {
    return integrated_[InterfaceDescription::VolumeFraction];
  }
  // interface states
  case InterfaceBlockBufferType::InterfaceStateVelocity: {
    return states_[InterfaceState::Velocity];
  }
  case InterfaceBlockBufferType::InterfaceStatePressurePositive: {
    return states_[InterfaceState::PressurePositive];
  }
  case InterfaceBlockBufferType::InterfaceStatePressureNegative: {
    return states_[InterfaceState::PressureNegative];
  }
  // interface parameters
#ifdef PERFORMANCE
  default: {
    return parameters_[InterfaceParameter::SurfaceTensionCoefficient];
  }
#else
  case InterfaceBlockBufferType::InterfaceParameterSurfaceTensionCoefficient: {
    return parameters_[InterfaceParameter::SurfaceTensionCoefficient];
  }
  default: {
    throw std::invalid_argument(
        "Requested buffer does not exist (impossible error");
  } break;
#endif
  }
}
