/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#ifndef FIELD_DETAILS_H
#define FIELD_DETAILS_H

#include <type_traits>
#include "utilities/helper_functions.h"
#include "enums/unit_type.h"

namespace FieldDetails {

   /**
   * @brief Unique identifier for all possible conservative equations in arbitrary order.
   * @note  Every member has to be added to the Equation enumeration as well.
   */
   enum class EquationPool {
      // mandatory equations
      Mass, Energy, MomentumX, MomentumY, MomentumZ,
      // optional equations
      // Example
   };

   /**
   * @brief Unique Identifier for all possible prime states in arbitrary order.
   * @note  Every member has to be added to the PrimeState enumeration as well.
   */
   enum class PrimeStatePool {
      // mandatory prime states
      Density, Pressure, VelocityX, VelocityY, VelocityZ,
      // optional prime states
      Temperature, // Example
   };

   /**
   * @brief Unique Identifier for all possible parameters in arbitrary order.
   * @note  Every member has to be added to the Parameters enumeration as well.
   */
   enum class ParameterPool {
      // viscosity parameters
      ShearViscosity,
      // conductivity parameters
      ThermalConductivity
   };

   /**
   * @brief Unique Identifier for all possible interface descriptions variables in arbitrary order.
   * @note  Every member has to be added to the InterfaceDescription enumeration as well.
   */
   enum class InterfaceDescriptionPool {
      Levelset, VolumeFraction
   };

   /**
   * @brief Unique Identifier for all possible interface states in arbitrary order.
   * @note  Every member has to be added to the InterfaceStates enumeration as well.
   */
   enum class InterfaceStatePool {
      // mandatory interface states
      Velocity,
      // optional interface states
      PressurePositive, PressureNegative
   };

   /**
   * @brief Unique Identifier for all possible interface parameters in arbitrary order.
   * @note  Every member has to be added to the InterfaceParameter enumeration as well.
   */
   enum class InterfaceParameterPool {
      SurfaceTensionCoefficient
   };

   /**
   * @brief Bundles relevant configuration data for active material fields (i.e. conservative equations, prime states, interface quantities).
   * @tparam DerivedActiveFieldsDefinition CRTP template parameter
   * @tparam FieldEnum Enumeration type of the underlying material field
   */
   template<typename DerivedActiveFieldsDefinition, typename FieldEnum>
   struct ActiveFieldsDefinition {
      struct FieldInfo {
         FieldEnum const Field;
         UnitType const Unit;
         std::string_view const InputName = "";
      };

      /**
      * @brief Value used to offset inactive from active material fields in the later generated enumeration.
      * @note IMPORTANT: Has to be larger than the number of active fields.
      */
      static constexpr unsigned int InactiveFieldOffset = 100;

      /**
      * @brief Gives the number of active fields.
      */
      static constexpr unsigned int Count = DerivedActiveFieldsDefinition::Definition.size();

      /**
      * @brief Gives the index of the given field. Should only be used to construct the resulting enumeration.
      * @param field Field whose index should be found.
      * @return Index of the field.
      */
      static constexpr unsigned int FieldIndex( FieldEnum const field ) {
         for( unsigned int i = 0; i < Count; ++i ) {
            if( DerivedActiveFieldsDefinition::Definition[i].Field == field ) return i;
         }
         return static_cast<unsigned int>( field ) + InactiveFieldOffset;
      }

      /**
      * @brief Gives the unit of the field at the given index.
      * @param index Index of the field.
      * @return UnitType representing the unit of the field.
      */
      static constexpr auto GetUnit( unsigned int const index ) {
         return DerivedActiveFieldsDefinition::Definition[index].Unit;
      }

      /**
      * @brief Gives the name of the field as it should be used in input reading.
      * @param index Index of the field.
      * @return Name used for input
      */
      static constexpr auto GetInputName( unsigned int const index ) {
         return DerivedActiveFieldsDefinition::Definition[index].InputName;
      }
   };

   /**
   * @brief Configuration of the active conservative equations.
   */
   struct ActiveEquations : public ActiveFieldsDefinition< ActiveEquations, EquationPool > {
      static constexpr auto Definition = MakeArray(
         //           conservative                  unit                    input_name*
         // mandatory equations (DO NOT CHANGE):
           FieldInfo{ EquationPool::Mass,           UnitType::Density,          "" }
         , FieldInfo{ EquationPool::Energy,         UnitType::Energy,           "" }
         , FieldInfo{ EquationPool::MomentumX,      UnitType::Momentum,         "" }
   #if DIMENSION != 1
         , FieldInfo{ EquationPool::MomentumY,      UnitType::Momentum,         "" }
   #endif
   #if DIMENSION == 3
         , FieldInfo{ EquationPool::MomentumZ,      UnitType::Momentum,         "" }
   #endif
         // optional equations:
      // , FieldInfo{ EquationPool::Example,        UnitType::Example,      "example" }
      );
      // notes: *  if input name is empty the field will not be read from input files. Spaces are not allowed
      static_assert( Definition.size() < InactiveFieldOffset, "Too many active equations! Please increase ActiveFieldsDefinition::InactiveFieldOffset." );
   };

   /**
   * @brief Configuration of the active prime states.
   */
   struct ActivePrimeStates : public ActiveFieldsDefinition< ActivePrimeStates, PrimeStatePool > {
      static constexpr auto Definition = MakeArray(
         //           prime state                      unit                       input_name
         // mandatory prime states (DO NOT CHANGE):
           FieldInfo{ PrimeStatePool::Density,         UnitType::Density,         "density" }
         , FieldInfo{ PrimeStatePool::Pressure,        UnitType::Pressure,        "pressure" }
         , FieldInfo{ PrimeStatePool::VelocityX,       UnitType::Velocity,        "velocityX" }
   #if DIMENSION != 1
         , FieldInfo{ PrimeStatePool::VelocityY,       UnitType::Velocity,        "velocityY" }
   #endif
   #if DIMENSION == 3
         , FieldInfo{ PrimeStatePool::VelocityZ,       UnitType::Velocity,        "velocityZ" }
   #endif
         // optional prime states:
         , FieldInfo{ PrimeStatePool::Temperature,     UnitType::Temperature,          "" }
      // , FieldInfo{ PrimeStatePool::Example,        UnitType::Example,      "example" }
      );
      // notes: *  if input name is empty the field will not be read from input files. Spaces are not allowed
      static_assert( Definition.size() < InactiveFieldOffset, "Too many active prime states! Please increase ActiveFieldsDefinition::InactiveFieldOffset." );
   };

   /**
   * @brief Configuration of the active parameters.
   */
   struct ActiveParameters : public ActiveFieldsDefinition< ActiveParameters, ParameterPool > {
      static constexpr auto Definition = MakeArray(
         //           parameter                                      unit                  input_name
           FieldInfo{ ParameterPool::ShearViscosity,         UnitType::Viscosity,            "" }
         , FieldInfo{ ParameterPool::ThermalConductivity,    UnitType::ThermalConductivity,  "" }

      // , FieldInfo{ ParameterPool::Example,                    UnitType::Example,      "example" }
      );
      // notes: *  if input name is empty the field will not be read from input files. Spaces are not allowed
      static_assert( Definition.size() < InactiveFieldOffset, "Too many active parameters! Please increase ActiveFieldsDefinition::InactiveFieldOffset." );
   };

   /**
   * @brief Configuration of the active interface states.
   */
   struct ActiveInterfaceDescriptions : public ActiveFieldsDefinition< ActiveInterfaceDescriptions, InterfaceDescriptionPool > {
      static constexpr auto Definition = MakeArray(
         //           interface description                               unit                  input_name
         // mandatory description (DO NOT CHANGE):
           FieldInfo{ InterfaceDescriptionPool::Levelset,                 UnitType::Unitless,     "phi" }
         , FieldInfo{ InterfaceDescriptionPool::VolumeFraction,      UnitType::Unitless,         "" }
      // , FieldInfo{ InterfaceDescriptionPool::Example,             UnitType::Example,      "example" }
      );
      // notes: *  if input name is empty the field will not be read from input files. Spaces are not allowed
      static_assert( Definition.size() < InactiveFieldOffset, "Too many active interface descriptions! Please increase ActiveFieldsDefinition::InactiveFieldOffset." );
   };

   /**
   * @brief Configuration of the active interface states.
   */
   struct ActiveInterfaceStates : public ActiveFieldsDefinition< ActiveInterfaceStates, InterfaceStatePool > {
      static constexpr auto Definition = MakeArray(
         //           interface state                             unit                  input_name
         // mandatory interface states (DO NOT CHANGE):
           FieldInfo{ InterfaceStatePool::Velocity,            UnitType::Velocity,         "" }
         // optional interface states:
         , FieldInfo{ InterfaceStatePool::PressurePositive,    UnitType::Pressure,         "" }
         , FieldInfo{ InterfaceStatePool::PressureNegative,    UnitType::Pressure,         "" }
      // , FieldInfo{ InterfaceStatePool::Example,             UnitType::Example,      "example" }
      );
      // notes: *  if input name is empty the field will not be read from input files. Spaces are not allowed
      static_assert( Definition.size() < InactiveFieldOffset, "Too many active interface states! Please increase ActiveFieldsDefinition::InactiveFieldOffset." );
   };

   /**
   * @brief Configuration of the active interface parameters.
   */
   struct ActiveInterfaceParameters : public ActiveFieldsDefinition< ActiveInterfaceParameters, InterfaceParameterPool > {
      static constexpr auto Definition = MakeArray(
         //           interface parameter                                                     unit                              input_name              component_name*
           FieldInfo{ InterfaceParameterPool::SurfaceTensionCoefficient,        UnitType::SurfaceTensionCoefficient,                "" }
      // , FieldInfo{ InterfaceParameterPool::Example,                                  UnitType::Example,                       "example" }
      );
      // notes: *  if input name is empty the field will not be read from input files. Spaces are not allowed
      static_assert( Definition.size() < InactiveFieldOffset, "Too many active interface parameters! Please increase ActiveFieldsDefinition::InactiveFieldOffset." );
   };
}

// NH: NEVER EVER change the underlying type.
/**
 * @brief Unique Identifier for the conservative equation. All active equations are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of EquationPool has to be added here as well.
 */
enum class Equation : unsigned int {
   // mandatory equations:
   Mass           = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::Mass ),
   Energy         = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::Energy ),
   MomentumX      = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::MomentumX ),
   MomentumY      = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::MomentumY ),
   MomentumZ      = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::MomentumZ ),
   // optional equations
   // Example     = FieldDetails::ActiveEquations::FieldIndex( FieldDetails::EquationPool::Example ),
};
/**
 * @brief Converts an equation identifier to a (C++11 standard compliant, i. e. positive) array index. "ETI = Equation to Index"
 * @param e The equation identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<Equation>::type ETI( Equation const e ) { return static_cast<typename std::underlying_type<Equation>::type>( e ); }

// NH: NEVER EVER change the underlying type.
/**
 * @brief Unique Identifier for the prime states. All active prime states are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of PrimeStatePool has to be added here as well.
 */
enum class PrimeState : unsigned int {
   // mandatory prime states
   Density        = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::Density ),
   Pressure       = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::Pressure ),
   VelocityX      = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::VelocityX ),
   VelocityY      = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::VelocityY ),
   VelocityZ      = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::VelocityZ ),
   // optional prime states
   Temperature    = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::Temperature ),
   // Example     = FieldDetails::ActivePrimeStates::FieldIndex( FieldDetails::PrimeStatePool::Example ),
};
/**
 * @brief Converts a prime state identifier to a (C++11 standard compliant, i. e. positive) array index. "PTI = Prime state to Index"
 * @param p The prime state identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<PrimeState>::type PTI( PrimeState const p ) { return static_cast<typename std::underlying_type<PrimeState>::type>( p ); }

// NH: NEVER EVER change the underlying type.
/**
 * @brief Unique Identifier for the parameters. All active parameters are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of ParameterPool has to be added here as well.
 */
enum class Parameter : unsigned int {
   // parameters
   ShearViscosity       = FieldDetails::ActiveParameters::FieldIndex( FieldDetails::ParameterPool::ShearViscosity ),
   ThermalConductivity  = FieldDetails::ActiveParameters::FieldIndex( FieldDetails::ParameterPool::ThermalConductivity )
   // Example     = FieldDetails::ActiveParameters::FieldIndex( FieldDetails::ParameterPool::Example ),
};

/**
 * @brief Converts a parameter identifier to a (C++11 standard compliant, i. e. positive) array index. "PTI = Parameter to Index"
 * @param p The parameter identifier.
 * @return Index to be used in Arrays.
 */
constexpr std::underlying_type<Parameter>::type PTI( Parameter const p ) { return static_cast<typename std::underlying_type<Parameter>::type>( p ); }

// NH: NEVER EVER change the underlying type.
/**
 * @brief Unique Identifier for the interface descriptions. All active interface descritpions are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of InterfaceDescriptionPool has to be added here as well.
 */
enum class InterfaceDescription : unsigned int {
   // mandatory interface descriptions
   Levelset        = FieldDetails::ActiveInterfaceDescriptions::FieldIndex( FieldDetails::InterfaceDescriptionPool::Levelset ),
   VolumeFraction  = FieldDetails::ActiveInterfaceDescriptions::FieldIndex( FieldDetails::InterfaceDescriptionPool::VolumeFraction ),
   // Example      = FieldDetails::ActiveInterfaceDescriptions::FieldIndex( FieldDetails::InterfaceDescriptionPool::Example ),
};

/**
 * @brief Converts a interface description identifier to a (C++11 standard compliant, i. e. positive) array index. "IDTI = Interface Descriptions To Index"
 * @param id The interface description identifier.
 * @return Index to be used in arrays.
 */
constexpr std::underlying_type<InterfaceDescription>::type IDTI( InterfaceDescription const id ) { return static_cast<typename std::underlying_type<InterfaceDescription>::type>( id ); }

// NH: NEVER EVER change the underlying type.
/**
 * @brief Unique Identifier for the interface states. All active interface states are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of InterfaceStatePool has to be added here as well.
 */
enum class InterfaceState : unsigned int {
   // mandatory interface states
   Velocity             = FieldDetails::ActiveInterfaceStates::FieldIndex( FieldDetails::InterfaceStatePool::Velocity ),
   // optional interface states
   PressurePositive     = FieldDetails::ActiveInterfaceStates::FieldIndex( FieldDetails::InterfaceStatePool::PressurePositive ),
   PressureNegative     = FieldDetails::ActiveInterfaceStates::FieldIndex( FieldDetails::InterfaceStatePool::PressureNegative ),
   // Example           = FieldDetails::ActiveInterfaceStates::FieldIndex( FieldDetails::InterfaceStatePool::Example ),
};

/**
 * @brief Converts a interface state identifier to a (C++11 standard compliant, i. e. positive) array index. "ISTI = Interface State To Index"
 * @param is The interface state identifier.
 * @return Index to be used in arrays.
 */
constexpr std::underlying_type<InterfaceState>::type ISTI( InterfaceState const is ) { return static_cast<typename std::underlying_type<InterfaceState>::type>( is ); }

// NH: NEVER EVER change the underlying type.
/**
 * @brief Unique Identifier for the interface states. All active interface states are guaranteed to be in consecutive order starting from 0.
 * @note  IMPORTANT: Every member of InterfaceStatePool has to be added here as well.
 */
enum class InterfaceParameter : unsigned int {
   // mandatory interface states
   SurfaceTensionCoefficient  = FieldDetails::ActiveInterfaceParameters::FieldIndex( FieldDetails::InterfaceParameterPool::SurfaceTensionCoefficient ),
   // Example           = FieldDetails::ActiveInterfaceParameters::FieldIndex( FieldDetails::InterfaceParameterPool::Example ),
};
/**
 * @brief Converts a interface parameter identifier to a (C++11 standard compliant, i. e. positive) array index. "ISTI = Interface Parameter To Index"
 * @param ip The interface parameter identifier.
 * @return Index to be used in arrays.
 */
constexpr std::underlying_type<InterfaceParameter>::type IPTI( InterfaceParameter const ip ) { return static_cast<typename std::underlying_type<InterfaceParameter>::type>( ip ); }

#endif // FIELD_DETAILS_H