//===------------------------- field_details.h ----------------------------===//
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
#ifndef FIELD_DETAILS_H
#define FIELD_DETAILS_H

#include "block_definitions/field_enums.h"
#include "block_definitions/field_setup.h"
#include "enums/unit_type.h"
#include "user_specifications/equation_settings.h"
#include "utilities/helper_functions.h"
#include <type_traits>

namespace FieldDetails {
/**
 * @brief Bundles relevant configuration data for active material fields (i.e.
 * conservative equations, prime states, interface quantities).
 * @tparam DerivedActiveFieldsDefinition CRTP template parameter.
 * @tparam FieldEnum Enumeration type of the underlying material field.
 */
template <typename DerivedActiveFieldsDefinition, typename FieldEnum>
struct ActiveFieldsDefinition {
  /**
   * @brief Value used to offset inactive from active material fields in the
   * later generated enumeration.
   * @note IMPORTANT: Has to be larger than the number of active fields.
   */
  static constexpr unsigned int inactive_field_offset = 100;

  /**
   * @brief Gives the number of active fields.
   */
  static constexpr unsigned int count =
      DerivedActiveFieldsDefinition::definition.size();

  /**
   * @brief Gives the index of the given field. Should only be used to construct
   * the resulting enumeration.
   * @param field Field whose index should be found.
   * @return Index of the field.
   */
  template <typename T = unsigned int>
  static constexpr T FieldIndex(FieldEnum const field) {
    for (unsigned int i = 0; i < count; ++i) {
      if (DerivedActiveFieldsDefinition::definition[i].Field == field)
        return static_cast<T>(i);
    }
    return static_cast<T>(static_cast<unsigned int>(field) +
                          inactive_field_offset);
  }

  /**
   * @brief Gives the unit of the field at the given index.
   * @param index Index of the field.
   * @return UnitType representing the unit of the field.
   */
  static constexpr auto GetUnit(unsigned int const index) {
    return DerivedActiveFieldsDefinition::definition[index].Unit;
  }

  /**
   * @brief Gives the name of the field as it should be used in input reading.
   * @param index Index of the field.
   * @return Name used for input.
   */
  static constexpr auto GetInputName(unsigned int const index) {
    return DerivedActiveFieldsDefinition::definition[index].Name;
  }
};

/**
 * @brief Configuration of the active conservative equations.
 */
struct ActiveEquations
    : public ActiveFieldsDefinition<ActiveEquations, EquationPool> {
  static constexpr auto definition =
      FieldSetup::EquationDefinition<active_equations>();
  static_assert(definition.size() < inactive_field_offset,
                "Too many active equations! Please increase "
                "ActiveFieldsDefinition::InactiveFieldOffset.");
};

/**
 * @brief Configuration of the active prime states.
 */
struct ActivePrimeStates
    : public ActiveFieldsDefinition<ActivePrimeStates, PrimeStatePool> {
  static constexpr auto definition =
      FieldSetup::PrimeStateDefinition<active_equations>();
  static_assert(definition.size() < inactive_field_offset,
                "Too many active prime states! Please increase "
                "ActiveFieldsDefinition::InactiveFieldOffset.");
};

/**
 * @brief Configuration of the active parameters.
 */
struct ActiveParameters
    : public ActiveFieldsDefinition<ActiveParameters, ParameterPool> {
  static constexpr auto definition =
      FieldSetup::ParameterDefinition<active_parameters>();
  static_assert(definition.size() < inactive_field_offset,
                "Too many active parameters! Please increase "
                "ActiveFieldsDefinition::inactive_field_offset.");
};

/**
 * @brief Configuration of the active interface description fields.
 */
struct ActiveInterfaceDescriptions
    : public ActiveFieldsDefinition<ActiveInterfaceDescriptions,
                                    InterfaceDescriptionPool> {
  static constexpr auto definition =
      FieldSetup::InterFaceDescriptionDefinition<active_interface_quantities>();
  static_assert(definition.size() < inactive_field_offset,
                "Too many active interface descriptions! Please increase "
                "ActiveFieldsDefinition::inactive_field_offset.");
};

/**
 * @brief Configuration of the active interface states.
 */
struct ActiveInterfaceStates
    : public ActiveFieldsDefinition<ActiveInterfaceStates, InterfaceStatePool> {
  static constexpr auto definition =
      FieldSetup::InterFaceStateDefinition<active_interface_quantities>();
  static_assert(definition.size() < inactive_field_offset,
                "Too many active interface quantities! Please increase "
                "ActiveFieldsDefinition::InactiveFieldOffset.");
};

/**
 * @brief Configuration of the active interface parameters.
 */
struct ActiveInterfaceParameters
    : public ActiveFieldsDefinition<ActiveInterfaceParameters,
                                    InterfaceParameterPool> {
  static constexpr auto definition =
      FieldSetup::InterFaceParameterDefinition<active_interface_quantities>();
  static_assert(definition.size() < inactive_field_offset,
                "Too many active interface parameters! Please increase "
                "ActiveFieldsDefinition::inactive_field_offset.");
};
} // namespace FieldDetails

/**
 * @brief Unique Identifier for the conservative equation. All active equations
 * are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of EquationPool has to be added here as well.
 */
enum class Equation : unsigned int {
  // mandatory equations:
  Mass = FieldDetails::ActiveEquations::FieldIndex(EquationPool::Mass),
  Energy = FieldDetails::ActiveEquations::FieldIndex(EquationPool::Energy),
  MomentumX =
      FieldDetails::ActiveEquations::FieldIndex(EquationPool::MomentumX),
  MomentumY =
      FieldDetails::ActiveEquations::FieldIndex(EquationPool::MomentumY),
  MomentumZ =
      FieldDetails::ActiveEquations::FieldIndex(EquationPool::MomentumZ),
  // optional equations
  Gamma = FieldDetails::ActiveEquations::FieldIndex(EquationPool::Gamma),
  Pi = FieldDetails::ActiveEquations::FieldIndex(EquationPool::Pi),
  // Example     = FieldDetails::ActiveEquations::FieldIndex(
  // EquationPool::Example ),
};
/**
 * @brief Converts an equation identifier to a (C++11 standard compliant, i. e.
 * positive) array index. "ETI = Equation To Index".
 * @param e The equation identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<Equation>::type ETI(Equation const e) {
  return static_cast<typename std::underlying_type<Equation>::type>(e);
}

/**
 * @brief Unique Identifier for the prime states. All active prime states are
 * guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of PrimeStatePool has to be added here as
 * well.
 */
enum class PrimeState : unsigned int {
  // mandatory prime states
  Density =
      FieldDetails::ActivePrimeStates::FieldIndex(PrimeStatePool::Density),
  Pressure =
      FieldDetails::ActivePrimeStates::FieldIndex(PrimeStatePool::Pressure),
  VelocityX =
      FieldDetails::ActivePrimeStates::FieldIndex(PrimeStatePool::VelocityX),
  VelocityY =
      FieldDetails::ActivePrimeStates::FieldIndex(PrimeStatePool::VelocityY),
  VelocityZ =
      FieldDetails::ActivePrimeStates::FieldIndex(PrimeStatePool::VelocityZ),
  // optional prime states
  Temperature =
      FieldDetails::ActivePrimeStates::FieldIndex(PrimeStatePool::Temperature),
  gamma = FieldDetails::ActivePrimeStates::FieldIndex(PrimeStatePool::gamma),
  pi = FieldDetails::ActivePrimeStates::FieldIndex(PrimeStatePool::pi),
  // Example     = FieldDetails::ActivePrimeStates::FieldIndex(
  // PrimeStatePool::Example ),
};
/**
 * @brief Converts a prime state identifier to a (C++11 standard compliant, i.
 * e. positive) array index. "PTI = Prime state To Index".
 * @param p The prime state identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<PrimeState>::type PTI(PrimeState const p) {
  return static_cast<typename std::underlying_type<PrimeState>::type>(p);
}

/**
 * @brief Unique Identifier for the parameters. All active parameters are
 * guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of ParameterPool has to be added here as well.
 */
enum class Parameter : unsigned int {
  // parameters
  ShearViscosity =
      FieldDetails::ActiveParameters::FieldIndex(ParameterPool::ShearViscosity),
  ThermalConductivity = FieldDetails::ActiveParameters::FieldIndex(
      ParameterPool::ThermalConductivity)
  // Example     = FieldDetails::ActiveParameters::FieldIndex(
  // ParameterPool::Example ),
};

/**
 * @brief Converts a parameter identifier to a (C++11 standard compliant, i. e.
 * positive) array index. "PTI = Parameter To Index".
 * @param p The parameter identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<Parameter>::type PTI(Parameter const p) {
  return static_cast<typename std::underlying_type<Parameter>::type>(p);
}

/**
 * @brief Unique Identifier for the interface descriptions. All active interface
 * descritpions are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of InterfaceDescriptionPool has to be added
 * here as well.
 */
enum class InterfaceDescription : unsigned int {
  // mandatory interface descriptions
  Levelset = FieldDetails::ActiveInterfaceDescriptions::FieldIndex(
      InterfaceDescriptionPool::Levelset),
  VolumeFraction = FieldDetails::ActiveInterfaceDescriptions::FieldIndex(
      InterfaceDescriptionPool::VolumeFraction),
  // Example      = FieldDetails::ActiveInterfaceDescriptions::FieldIndex(
  // InterfaceDescriptionPool::Example ),
};

/**
 * @brief Converts a interface description identifier to a (C++11 standard
 * compliant, i. e. positive) array index. "IDTI = Interface Descriptions To
 * Index".
 * @param id The interface description identifier.
 * @return Index to be used in arrays.
 */
constexpr std::underlying_type<InterfaceDescription>::type
IDTI(InterfaceDescription const id) {
  return static_cast<typename std::underlying_type<InterfaceDescription>::type>(
      id);
}

/**
 * @brief Unique Identifier for the interface states. All active interface
 * states are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of InterfaceStatePool has to be added here as
 * well.
 */
enum class InterfaceState : unsigned int {
  // mandatory interface states
  Velocity = FieldDetails::ActiveInterfaceStates::FieldIndex(
      InterfaceStatePool::Velocity),
  // optional interface states
  PressurePositive = FieldDetails::ActiveInterfaceStates::FieldIndex(
      InterfaceStatePool::PressurePositive),
  PressureNegative = FieldDetails::ActiveInterfaceStates::FieldIndex(
      InterfaceStatePool::PressureNegative),
  // Example           = FieldDetails::ActiveInterfaceStates::FieldIndex(
  // InterfaceStatePool::Example ),
};

/**
 * @brief Converts a interface state identifier to a (C++11 standard compliant,
 * i. e. positive) array index. "ISTI = Interface State To Index".
 * @param is The interface state identifier.
 * @return Index to be used in arrays.
 */
constexpr std::underlying_type<InterfaceState>::type
ISTI(InterfaceState const is) {
  return static_cast<typename std::underlying_type<InterfaceState>::type>(is);
}

/**
 * @brief Unique Identifier for the interface states. All active interface
 * states are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of InterfaceStatePool has to be added here as
 * well.
 */
enum class InterfaceParameter : unsigned int {
  // mandatory interface states
  SurfaceTensionCoefficient =
      FieldDetails::ActiveInterfaceParameters::FieldIndex(
          InterfaceParameterPool::SurfaceTensionCoefficient),
  // Example           = FieldDetails::ActiveInterfaceParameters::FieldIndex(
  // InterfaceParameterPool::Example ),
};
/**
 * @brief Converts a interface parameter identifier to a (C++11 standard
 * compliant, i. e. positive) array index. "ISTI = Interface Parameter To
 * Index".
 * @param ip The interface parameter identifier.
 * @return Index to be used in arrays.
 */
constexpr std::underlying_type<InterfaceParameter>::type
IPTI(InterfaceParameter const ip) {
  return static_cast<typename std::underlying_type<InterfaceParameter>::type>(
      ip);
}

#endif // FIELD_DETAILS_H
