//===------------------------ field_materials.h ---------------------------===//
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
#ifndef FIELD_MATERIAL_DEFINITONS_H
#define FIELD_MATERIAL_DEFINITONS_H

#include "block_definitions/field_details.h"
#include "user_specifications/compile_time_constants.h"
#include "utilities/container_operations.h"
#include "utilities/helper_functions.h"
#include <algorithm>
#include <array>
#include <tuple>
#include <type_traits>

/**
 * @brief Unique identifier for the conservative buffer type, i.e. the average,
 * right-hand side, or initial buffer.
 */
enum class ConservativeBufferType { Average, RightHandSide, Initial };

/**
 * @brief Unique Identifier for the material Field type, i.e. a conservative, a
 * prime-state or a parameter field.
 */
enum class MaterialFieldType { Conservatives, PrimeStates, Parameters };

class MaterialFieldsDefinitions {

  // get arrays of the consecutively ordered active fields
  static constexpr auto active_equations_ = IndexSequenceToEnumArray<Equation>(
      std::make_index_sequence<FieldDetails::ActiveEquations::count>{});
  static constexpr auto active_prime_states_ =
      IndexSequenceToEnumArray<PrimeState>(
          std::make_index_sequence<FieldDetails::ActivePrimeStates::count>{});
  static constexpr auto active_parameters_ =
      IndexSequenceToEnumArray<Parameter>(
          std::make_index_sequence<FieldDetails::ActiveParameters::count>{});

  static constexpr std::array<Equation, DTI(CC::DIM())> const active_momenta_ =
      {Equation::MomentumX
#if DIMENSION != 1
       ,
       Equation::MomentumY
#endif
#if DIMENSION == 3
       ,
       Equation::MomentumZ
#endif
  };

  static constexpr std::array<PrimeState, DTI(CC::DIM())> const
      active_velocities_ = {PrimeState::VelocityX
#if DIMENSION != 1
                            ,
                            PrimeState::VelocityY
#endif
#if DIMENSION == 3
                            ,
                            PrimeState::VelocityZ
#endif
  };

  static constexpr auto wavelet_analysis_equations_ =
      ContainerOperations::ArrayOfElementWiseFunctionApplication(
          FieldSetup::WaveletEquations<active_equations>(),
          FieldDetails::ActiveEquations::FieldIndex<Equation>);

public:
  MaterialFieldsDefinitions() = delete;
  ~MaterialFieldsDefinitions() = default;
  MaterialFieldsDefinitions(MaterialFieldsDefinitions const &) = delete;
  MaterialFieldsDefinitions &
  operator=(MaterialFieldsDefinitions const &) = delete;
  MaterialFieldsDefinitions(MaterialFieldsDefinitions &&) = delete;
  MaterialFieldsDefinitions &operator=(MaterialFieldsDefinitions &&) = delete;

  /**
   * @brief Gives whether the given field name and index is active.
   * @param field_type Identifier for the material field type (Conservatives,
   * PrimeStates or Parameters).
   * @param field_index Index of the given field.
   * @return True if active, otherwise False.
   */
  static constexpr bool IsFieldActive(MaterialFieldType const field_type,
                                      unsigned int const field_index) {
    switch (field_type) {
    case MaterialFieldType::Conservatives: {
      return field_index < FieldDetails::ActiveEquations::inactive_field_offset;
    }
    case MaterialFieldType::Parameters: {
      return field_index <
             FieldDetails::ActiveParameters::inactive_field_offset;
    }
    default: { // MaterialFieldType::PrimeStates
      return field_index <
             FieldDetails::ActivePrimeStates::inactive_field_offset;
    }
    }
  }

  /**
   * @brief Gives whether the given conservative equation is active.
   * @param eq Equation identifier.
   * @return True if active, false otherwise.
   */
  static constexpr bool IsEquationActive(Equation const eq) {
    return static_cast<unsigned int>(eq) <
           FieldDetails::ActiveEquations::inactive_field_offset;
  }

  /**
   * @brief Gives whether the given prime state is active.
   * @param ps PrimeState identifier.
   * @return True if active, false otherwise.
   */
  static constexpr bool IsPrimeStateActive(PrimeState const ps) {
    return static_cast<unsigned int>(ps) <
           FieldDetails::ActivePrimeStates::inactive_field_offset;
  }

  /**
   * @brief Gives whether the given parameter is active.
   * @param pa Parameter identifier.
   * @return True if active, false otherwise.
   */
  static constexpr bool IsParameterActive(Parameter const pa) {
    return static_cast<unsigned int>(pa) <
           FieldDetails::ActiveParameters::inactive_field_offset;
  }

  /**
   * @brief Gives the name used for the reading of input files for the given
   * field and index .
   * @param field_type Identifier for the material field type (Conservatives,
   * PrimeStates or Parameters).
   * @param field_index Index of the given field.
   * @return Input name for the given field and index (empty if not specified).
   */
  static constexpr auto InputName(MaterialFieldType const field_type,
                                  unsigned int const field_index) {
    switch (field_type) {
    case MaterialFieldType::Conservatives: {
      return FieldDetails::ActiveEquations::GetInputName(field_index);
    }
    case MaterialFieldType::Parameters: {
      return FieldDetails::ActiveParameters::GetInputName(field_index);
    }
    default: { // MaterialFieldType::PrimeStates
      return FieldDetails::ActivePrimeStates::GetInputName(field_index);
    }
    }
  }

  /**
   * @brief Gives the name used in the output for the given equation.
   * @param eq Equation identifier.
   * @return Input name for the given equation (empty if not specified).
   */
  static constexpr auto InputName(Equation const eq) {
    return FieldDetails::ActiveEquations::GetInputName(ETI(eq));
  }

  /**
   * @brief Gives the name used in the output for the given prime state.
   * @param ps PrimeState identifier.
   * @return Input name for the given primestate (empty if not specified).
   */
  static constexpr auto InputName(PrimeState const ps) {
    return FieldDetails::ActivePrimeStates::GetInputName(PTI(ps));
  }

  /**
   * @brief Gives the name used in the output for the given parameter.
   * @param pa Parameter identifier.
   * @return Input name for the given parameter (empty if not specified).
   */
  static constexpr auto InputName(Parameter const pa) {
    return FieldDetails::ActiveParameters::GetInputName(PTI(pa));
  }

  /**
   * @brief Gives the dimension/unit for the given field and index.
   * @param field_type Identifier for the material field type (Conservatives,
   * PrimeStates or Parameters).
   * @param field_index Index of the given field.
   * @return Unit for the given field and index.
   */
  static constexpr auto FieldUnit(MaterialFieldType const field_type,
                                  unsigned int const field_index) {
    switch (field_type) {
    case MaterialFieldType::Conservatives: {
      return FieldDetails::ActiveEquations::GetUnit(field_index);
    }
    case MaterialFieldType::Parameters: {
      return FieldDetails::ActiveParameters::GetUnit(field_index);
    }
    default: { // MaterialFieldType::PrimeStates
      return FieldDetails::ActivePrimeStates::GetUnit(field_index);
    }
    }
  }

  /**
   * @brief Gives the dimension/unit for the given equation.
   * @param eq Equation identifier.
   * @return Unit for the given equation.
   */
  static constexpr auto FieldUnit(Equation const eq) {
    return FieldDetails::ActiveEquations::GetUnit(ETI(eq));
  }

  /**
   * @brief Gives the dimension/unit for the given prime state.
   * @param ps PrimeState identifier.
   * @return Unit for the given prime state.
   */
  static constexpr auto FieldUnit(PrimeState const ps) {
    return FieldDetails::ActivePrimeStates::GetUnit(PTI(ps));
  }

  /**
   * @brief Gives the dimension/unit for the given parameter.
   * @param pa Parameter identifier.
   * @return Unit for the given parameter.
   */
  static constexpr auto FieldUnit(Parameter const pa) {
    return FieldDetails::ActiveParameters::GetUnit(PTI(pa));
  }

  /**
   * @brief Gives the number of equations considered in the simulation, i.e.
   * Euler Equations.
   * @return Number of Equations = 5 (mass, energy, X-,Y-,Z-momentum) for 3D.
   * @return Number of Equations = 4 (mass, energy, X-,Y-momentum) for 2D.
   * @return Number of Equations = 3 (mass, energy, X-momentum) for 1D.
   *
   * @note Depending on the configuration of active conservatives, the number
   * can change.
   */
  static constexpr unsigned int ANOE() {
    return FieldDetails::ActiveEquations::count;
  }

  /**
   * @brief Gives the number of prime states considered in the simulation.
   * @return Number of prime states = 6 (density, pressure, X-, Y-, Z-velocity)
   * for 3D.
   * @return Number of prime states = 5 (density, pressure, X-, Y-velocity) for
   * 2D.
   * @return Number of prime states = 4 (density, pressure, X-velocity) for 1D.
   *
   * @note Depending on the configuration of active prime states, the number can
   * change (e.g., temperature activated).
   */
  static constexpr unsigned int ANOP() {
    return FieldDetails::ActivePrimeStates::count;
  }

  /**
   * @brief Gives the number of parameters considered in the simulation.
   * @return Number of active parameters
   *
   * @note Depending on the configuration of active parameters, the number can
   * change.
   */
  static constexpr unsigned int ANOPA() {
    return FieldDetails::ActiveParameters::count;
  }

  /**
   * @brief Gives the number of active fields for the given field type, i.e.
   * conservatives, prime states or parameters.
   * @param field_type The material-field type.
   * @return Number of active fields.
   */
  static constexpr unsigned int ANOF(MaterialFieldType const field_type) {
    switch (field_type) {
    case MaterialFieldType::Conservatives: {
      return ANOE();
    }
    case MaterialFieldType::Parameters: {
      return ANOPA();
    }
    default: { // MaterialFieldType::PrimeStates
      return ANOP();
    }
    }
  }

  /**
   * @brief Gives the set of equations which are worked on in this simulation.
   * I. e. varies with dimension of the simulation. "ASOE = Active Set of
   * Equations".
   * @return Set of equations. E.g. Rho, Energy, X-Momentum for a 1D pure
   * material simulation.
   */
  static constexpr auto ASOE() { return active_equations_; }

  /**
   * @brief Gives the active momentum equations, i.e. varies with dimension of
   * the simulation. "AME = Active Momentum Equations".
   * @return Set of momentum equations.
   */
  static constexpr auto AME() { return active_momenta_; }

  /**
   * @brief Gives the set of prime states which are worked on in this
   * simulation. I. e. varies with dimension of the simulation. "ASOP = Active
   * Set of Prime states".
   * @return Set of prime states. E.g. Rho, Pressure, X-Velocity for a 1D pure
   * material simulation.
   */
  static constexpr auto ASOP() { return active_prime_states_; }

  /**
   * @brief Gives the active velocity prime states, i.e. varies with dimension
   * of the simulation. "AV = Active Velocities".
   * @return Set of velocity prime states.
   */
  static constexpr auto AV() { return active_velocities_; }

  /**
   * @brief Gives the equations considered for the coarsening/refinement
   * decision. "EWA = Equations for Wavelet-Analysis".
   * @return List of equations.
   */
  static constexpr auto EWA() { return wavelet_analysis_equations_; }

  /**
   * @brief Gives the set of parameters which are used in this simulation "ASOPA
   * = Active Set of PArameters".
   * @return Set of parameters.
   */
  static constexpr auto ASOPA() { return active_parameters_; }

  /**
   * @brief Helper function to get the nth conservative equation that is NOT a
   * momentum equation.
   * @param n
   * @return Identifier of the nth non-momentum equation.
   */
  static constexpr Equation NthNonMomentumEquation(unsigned int const n) {
    unsigned int skipped = 0;
    for (unsigned int i = 0; i <= n;) {
      if (active_equations_[i + skipped] == Equation::MomentumX ||
          active_equations_[i + skipped] == Equation::MomentumY ||
          active_equations_[i + skipped] == Equation::MomentumZ) {

        ++skipped;
      } else {
        ++i;
      }
    }
    return active_equations_[n + skipped];
  }
};

using MF = MaterialFieldsDefinitions;

static_assert(std::make_pair(MF::AME()[0], MF::AV()[0]) ==
                  std::make_pair(Equation::MomentumX, PrimeState::VelocityX),
              "MF::AME()[0] and MF::AV()[0] have to consistently return "
              "MomentumX and VelocityX");
#if DIMENSION != 1
static_assert(std::make_pair(MF::AME()[1], MF::AV()[1]) ==
                  std::make_pair(Equation::MomentumY, PrimeState::VelocityY),
              "MF::AME()[1] and MF::AV()[1] have to consistently return "
              "MomentumY and VelocityY");
#endif
#if DIMENSION == 3
static_assert(std::make_pair(MF::AME()[2], MF::AV()[2]) ==
                  std::make_pair(Equation::MomentumZ, PrimeState::VelocityZ),
              "MF::AME()[2] and MF::AV()[2] have to consistently return "
              "MomentumZ and VelocityZ");
#endif

#endif // FIELDS_MATERIAL_DEFINITONS_H
