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
#ifndef FIELD_INTERFACE_DEFINITIONS_H
#define FIELD_INTERFACE_DEFINITIONS_H

#include <algorithm>
#include <array>
#include <tuple>
#include <type_traits>
#include "user_specifications/compile_time_constants.h"
#include "block_definitions/field_details.h"

/**
 * @brief Unique identifier for the interface description buffer type, i.e. the base, right-hand side, reinitialized or initial buffer.
 */
enum class InterfaceDescriptionBufferType { Base,
                                            RightHandSide,
                                            Reinitialized,
                                            Initial,
                                            Integrated };

/**
 * @brief Unique Identifier for the interface field type, i.e. a interface description, a interface state or a interface parameter field.
 */
enum class InterfaceFieldType { Description,
                                States,
                                Parameters };

class InterfaceFieldsDefinitions {

   // get arrays of the consecutively ordered active fields
   static constexpr auto active_description_ = IndexSequenceToEnumArray<InterfaceDescription>( std::make_index_sequence<FieldDetails::ActiveInterfaceDescriptions::count>{} );
   static constexpr auto active_states_      = IndexSequenceToEnumArray<InterfaceState>( std::make_index_sequence<FieldDetails::ActiveInterfaceStates::count>{} );
   static constexpr auto active_parameters_  = IndexSequenceToEnumArray<InterfaceParameter>( std::make_index_sequence<FieldDetails::ActiveInterfaceParameters::count>{} );
   // definition of interface states and parameters (the indices) which should be extended
   static constexpr std::array<InterfaceState, 1> states_to_extend_         = { InterfaceState::Velocity };
   static constexpr std::array<InterfaceParameter, 1> parameters_to_extend_ = { InterfaceParameter::SurfaceTensionCoefficient };

public:
   InterfaceFieldsDefinitions()                                    = delete;
   ~InterfaceFieldsDefinitions()                                   = default;
   InterfaceFieldsDefinitions( InterfaceFieldsDefinitions const& ) = delete;
   InterfaceFieldsDefinitions& operator=( InterfaceFieldsDefinitions const& ) = delete;
   InterfaceFieldsDefinitions( InterfaceFieldsDefinitions&& )                 = delete;
   InterfaceFieldsDefinitions& operator=( InterfaceFieldsDefinitions&& ) = delete;

   /**
    * @brief Gives whether the given field name and index is active.
    * @param field_type Identifier for the interface field type (Description, States or Parameters).
    * @param field_index Index of the given field.
    * @return True if active, otherwise False.
    */
   static constexpr bool IsFieldActive( InterfaceFieldType const field_type, unsigned int const field_index ) {
      switch( field_type ) {
         case InterfaceFieldType::Description: {
            return field_index < FieldDetails::ActiveInterfaceDescriptions::inactive_field_offset;
         }
         case InterfaceFieldType::Parameters: {
            return field_index < FieldDetails::ActiveInterfaceParameters::inactive_field_offset;
         }
         default: {// InterfaceFieldType::States
            return field_index < FieldDetails::ActiveInterfaceStates::inactive_field_offset;
         }
      }
   }

   /**
    * @brief Gives whether the given description is active.
    * @param d Description identifier.
    * @return True if active, false otherwise.
    */
   static constexpr bool IsDescriptionActive( InterfaceDescription const id ) {
      return static_cast<unsigned int>( id ) < FieldDetails::ActiveInterfaceDescriptions::inactive_field_offset;
   }

   /**
    * @brief Gives whether the given interface state is active.
    * @param is State identifier.
    * @return True if active, false otherwise.
    */
   static constexpr bool IsStateActive( InterfaceState const is ) {
      return static_cast<unsigned int>( is ) < FieldDetails::ActiveInterfaceStates::inactive_field_offset;
   }

   /**
    * @brief Gives whether the given parameter is active.
    * @param pa parameter identifier.
    * @return True if active, false otherwise.
    */
   static constexpr bool IsParameterActive( InterfaceParameter const pa ) {
      return static_cast<unsigned int>( pa ) < FieldDetails::ActiveInterfaceParameters::inactive_field_offset;
   }

   /**
    * @brief Gives the name used for reading input files for the given field and index.
    * @param field_type Interface field type for which the output name should be returned.
    * @param field_index field index for given type for which the output name should be returned.
    * @return Input name for field and index (empty if not specified).
    */
   static constexpr auto InputName( InterfaceFieldType const field_type, unsigned int const field_index ) {
      switch( field_type ) {
         case InterfaceFieldType::Description: {
            return FieldDetails::ActiveInterfaceDescriptions::GetInputName( field_index );
         }
         case InterfaceFieldType::Parameters: {
            return FieldDetails::ActiveInterfaceParameters::GetInputName( field_index );
         }
         default: {// InterfaceFieldType::States
            return FieldDetails::ActiveInterfaceStates::GetInputName( field_index );
         }
      }
   }

   /**
    * @brief Gives the name used in the input for the given field.
    * @param id InterfaceDescription type for which the input name should be returned.
    * @return Input name for interface description field (empty if not specified)
    */
   static constexpr auto InputName( InterfaceDescription const id ) {
      return FieldDetails::ActiveInterfaceDescriptions::GetInputName( IDTI( id ) );
   }

   /**
    * @brief Gives the name used in the input for the given field.
    * @param is Interface state type for which the input name should be returned
    * @return Input name for interface state field (empty if not specified)
    */
   static constexpr auto InputName( InterfaceState const is ) {
      return FieldDetails::ActiveInterfaceStates::GetInputName( ISTI( is ) );
   }

   /**
    * @brief Gives the name used in the input for the given field.
    * @param ip Interface parameter type for which the input name should be returned
    * @return Input name for interface parameter field (empty if not specified)
    */
   static constexpr auto InputName( InterfaceParameter const ip ) {
      return FieldDetails::ActiveInterfaceParameters::GetInputName( IPTI( ip ) );
   }

   /**
    * @brief Gives the dimension/unit of the given field and index
    * @param field_type Interface field type for which the unit/dimension should be returned
    * @param field_index field index for given type for which the unit/dimension should be returned
    * @return unit/dimension for field and index
    */
   static constexpr auto FieldUnit( InterfaceFieldType const field_type, unsigned int const field_index ) {
      switch( field_type ) {
         case InterfaceFieldType::Description: {
            return FieldDetails::ActiveInterfaceDescriptions::GetUnit( field_index );
         }
         case InterfaceFieldType::Parameters: {
            return FieldDetails::ActiveInterfaceParameters::GetUnit( field_index );
         }
         default: {// InterfaceFieldType::States
            return FieldDetails::ActiveInterfaceStates::GetUnit( field_index );
         }
      }
   }

   /**
    * @brief Gives the dimension/unit of the given field.
    * @param id InterfaceDescription type for which the unit should be returned.
    * @return Unit of desired field.
    */
   static constexpr auto FieldUnit( InterfaceDescription const id ) {
      return FieldDetails::ActiveInterfaceDescriptions::GetUnit( IDTI( id ) );
   }

   /**
    * @brief Gives the dimension/unit of the given field.
    * @param is Interface State type for which the unit should be returned.
    * @return Unit of desired field.
    */
   static constexpr auto FieldUnit( InterfaceState const is ) {
      return FieldDetails::ActiveInterfaceStates::GetUnit( ISTI( is ) );
   }

   /**
    * @brief Gives the dimension/unit of the given field.
    * @param ip Interface Parameter type for which the unit should be returned.
    * @return Unit of desired field.
    */
   static constexpr auto FieldUnit( InterfaceParameter const ip ) {
      return FieldDetails::ActiveInterfaceParameters::GetUnit( IPTI( ip ) );
   }

   /**
    * @brief Gives the number of InterfaceDescription considered in the simulation.
    * @return Number of active InterfaceDescription.
    *
    * @note Depending on the configuration of active inteface description, the number can change.
    */
   static constexpr unsigned int ANOD() { return FieldDetails::ActiveInterfaceDescriptions::count; }

   /**
    * @brief Gives the number of interface states considered in the simulation.
    * @return Number of active interface states.
    *
    * @note Depending on the configuration of active prime states, the number can change (e.g., pressure activated).
    */
   static constexpr unsigned int ANOS() { return FieldDetails::ActiveInterfaceStates::count; }

   /**
    * @brief Gives the number of interface parameters considered in the simulation.
    * @return Number of active interface parameters.
    *
    * @note Depending on the configuration of active prime states, the number can change.
    */
   static constexpr unsigned int ANOPA() { return FieldDetails::ActiveInterfaceParameters::count; }

   /**
    * @brief Gives the number of active fields for the given field type, i.e. description, states or parameters.
    * @param field_type The interface field type.
    * @return Number of active fields.
    */
   static constexpr auto ANOF( InterfaceFieldType const field_type ) {
      switch( field_type ) {
         case InterfaceFieldType::Description: {
            return ANOD();
         }
         case InterfaceFieldType::Parameters: {
            return ANOPA();
         }
         default: {// InterfaceFieldType::States
            return ANOS();
         }
      }
   }

   /**
    * @brief Gives the set of Description which are worked on in this simulation. "ASOD = Active Set of Description".
    * @return Set of Description. E.g. levelset and volume fraction.
    */
   static constexpr auto ASOD() { return active_description_; }

   /**
    * @brief Gives the set of interface states which are worked on in this simulation. "ASOS = Active Set of States".
    * @return Set of interface states.
    */
   static constexpr auto ASOS() { return active_states_; }

   /**
    * @brief Gives the set of interface parameters which are worked on in this simulation. "ASOPA = Active Set of PArameters".
    * @return Set of interface parameters.
    */
   static constexpr auto ASOPA() { return active_parameters_; }

   /**
    * @brief Gives the field to the appropriate field_index for the Description. "ITID = Index to Interface Description".
    * @param description_index Index in the interface description field for which the type should be returned.
    * @return desired interface description field.
    */
   static constexpr auto ITID( unsigned int const description_index ) {
      return active_description_[description_index];
   }

   /**
    * @brief Gives the field to the appropriate field_index for the interface states. "ITIS = Index to Interface State".
    * @param state_index Index in the state field for which the type should be returned.
    * @return desired state field.
    */
   static constexpr auto ITIS( unsigned int const state_index ) {
      return active_states_[state_index];
   }

   /**
    * @brief Gives the field to the appropriate field_index for the interface parameters. "ITIP = Index To Interface Parameter".
    * @param parameter_index Index in the parameter field for which the type should be returned.
    * @return desired parameter field.
    */
   static constexpr auto ITIP( unsigned int const parameter_index ) {
      return active_parameters_[parameter_index];
   }

   /**
    * @brief Gives the interface states which have to be extended at the interface. "ISTE = Interface States To Extend".
    * @return Interface states which have to be extended at the interface.
    */
   static constexpr auto ISTE() { return states_to_extend_; }

   /**
    * @brief Gives the number of interface states which have to be extended at the interface. "NOSTE = Number Of States To Extend".
    * @return Number of interface states which have to be extended at the interface.
    */
   static constexpr unsigned int NOSTE() { return states_to_extend_.size(); }

   /**
    * @brief Gives the number of interface parameters which have to be extended at the interface. "NOPATE = Number Of Parameters To Extend".
    * @return Number of interface parameters which have to be extended at the interface.
    */
   static constexpr unsigned int NOPATE() { return parameters_to_extend_.size(); }

   /**
    * @brief Gives the number of interface fields to be extended (for states or parameters). "NOFTE = Number of fields to extend".
    * @param field_type The interface field type.
    * @return Number of fields to be extended.
    */
   static constexpr auto NOFTE( InterfaceFieldType const field_type ) {
      switch( field_type ) {
         case InterfaceFieldType::Parameters: {
            return NOPATE();
         }
#ifndef PERFORMANCE
         case InterfaceFieldType::States: {
            return NOSTE();
         }
         default: {
            throw std::logic_error( "Interface field type not known!" );
         }
#else
         default: {// InterfaceFieldType::States
            return NOSTE();
         }
#endif
      }
   }

   /**
    * @brief Gives the appropriate index to the given interface field to be extended (for states or parameters). "FITE = Field Index to extend".
    * @param field_type The interface field type.
    * @param field_index Index of the field.
    * @return Appropriate index of the desired field in the active field definitions.
    *
    * @note Since no sanity check is done if the index is valid this function should ALWAYS be called in conjunction with the NOFTE() function
    *       that gives the appropriate number of fields that should be extended
    */
   static constexpr auto FITE( InterfaceFieldType const field_type, unsigned int const field_index ) {
      switch( field_type ) {
         case InterfaceFieldType::Parameters: {
            return IPTI( parameters_to_extend_[field_index] );
         }
#ifndef PERFORMANCE
         case InterfaceFieldType::States: {
            return ISTI( states_to_extend_[field_index] );
         }
         default: {
            throw std::logic_error( "Interface field type not known!" );
         }
#else
         default: {// InterfaceFieldType::States
            return ISTI( states_to_extend_[field_index] );
         }
#endif
      }
   }
};

using IF = InterfaceFieldsDefinitions;

#endif// FIELD_INTERFACE_DEFINITIONS_H
