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
#ifndef FIELD_MATERIAL_DEFINITONS_H
#define FIELD_MATERIAL_DEFINITONS_H

#include <algorithm>
#include <array>
#include <tuple>
#include <type_traits>
#include "user_specifications/compile_time_constants.h"
#include "field_details.h"

/**
 * @brief Unique identifier for the conservative buffer type, i.e. the average, right-hand side, or initial buffer.
 */
enum class ConservativeBufferType : unsigned short { Average = 0, RightHandSide = 1, Initial = 2 };

/**
 * @brief Converts an conservative buffer identifier to a (C++11 standard compliant, i. e. positive) array index. "CTTI = Conservative buffer to index"
 * @param cb The conservative buffer identifier
 * @return The index.
 */
constexpr std::underlying_type<ConservativeBufferType>::type CBTI( ConservativeBufferType const cb ) { return static_cast<typename std::underlying_type<ConservativeBufferType>::type>( cb ); }

/**
 * @brief Unique Identifier for the material Field type, i.e. a conservative, a prime-state or a parameter field.
 */
enum class MaterialFieldType : unsigned short { Conservatives = 0, PrimeStates = 1, Parameters = 2 };

/**
 * @brief Converts an material field identifier to a (C++11 standard compliant, i. e. positive) array index. "MFTI = Material Field to Index"
 * @param fft The material field type identifier
 * @return The index.
 */
constexpr std::underlying_type<MaterialFieldType>::type MFTI( MaterialFieldType const fft ) { return static_cast<typename std::underlying_type<MaterialFieldType>::type>( fft ); }


class MaterialFieldsDefinitions {

   // get arrays of the consecutively ordered active fields
   static constexpr auto active_equations_    = IndexSequenceToEnumArray<Equation>( std::make_index_sequence<FieldDetails::ActiveEquations::Count>{} );
   static constexpr auto active_prime_states_ = IndexSequenceToEnumArray<PrimeState>( std::make_index_sequence<FieldDetails::ActivePrimeStates::Count>{} );
   static constexpr auto active_parameters_   = IndexSequenceToEnumArray<Parameter>( std::make_index_sequence<FieldDetails::ActiveParameters::Count>{} );

   static constexpr const std::array<Equation,DTI( CC::DIM() )> active_momenta_ = { Equation::MomentumX
#if DIMENSION!=1
      ,Equation::MomentumY
#endif
#if DIMENSION==3
      ,Equation::MomentumZ
#endif
   };

   static constexpr const std::array<PrimeState,DTI( CC::DIM() )> active_velocities_ = { PrimeState::VelocityX
#if DIMENSION!=1
      ,PrimeState::VelocityY
#endif
#if DIMENSION==3
      ,PrimeState::VelocityZ
#endif
   };
   static constexpr std::array<Equation, 2> wavelet_analysis_equations_ { Equation::Mass, Equation::Energy };

public:
   MaterialFieldsDefinitions() = delete;
   ~MaterialFieldsDefinitions() = default;
   MaterialFieldsDefinitions( MaterialFieldsDefinitions const& ) = delete;
   MaterialFieldsDefinitions& operator=( MaterialFieldsDefinitions const& ) = delete;
   MaterialFieldsDefinitions( MaterialFieldsDefinitions&& ) = delete;
   MaterialFieldsDefinitions& operator=( MaterialFieldsDefinitions&& ) = delete;

   /**
    * @brief Gives whether the given field name and index is active
    * @param field_type Identifier for the material field type (Conservatives, PrimeStates or Parameters)
    * @param field_index Index of the given field 
    * @return True if active, otherwise False
    */ 
   static constexpr bool IsFieldActive( MaterialFieldType const field_type, unsigned int const field_index ) {
      switch( field_type ) {
         case MaterialFieldType::Conservatives : {
            return field_index < FieldDetails::ActiveEquations::InactiveFieldOffset;
         }
         case MaterialFieldType::Parameters : {
            return field_index < FieldDetails::ActiveParameters::InactiveFieldOffset;
         }
         default : { // MaterialFieldType::PrimeStates
            return field_index < FieldDetails::ActivePrimeStates::InactiveFieldOffset;
         }
      }
   }

   /**
    * @brief Gives whether the given conservative equation is active.
    * @param eq Equation identifier.
    * @return True if active, false otherwise.
    */
   static constexpr bool IsEquationActive( Equation const eq ) {
      return static_cast<unsigned int>( eq ) < FieldDetails::ActiveEquations::InactiveFieldOffset;
   }

   /**
    * @brief Gives whether the given prime state is active.
    * @param ps PrimeState identifier.
    * @return True if active, false otherwise.
    */
   static constexpr bool IsPrimeStateActive( PrimeState const ps ) {
      return static_cast<unsigned int>( ps ) < FieldDetails::ActivePrimeStates::InactiveFieldOffset;
   }

   /**
    * @brief Gives whether the given parameter is active.
    * @param pa Parameter identifier.
    * @return True if active, false otherwise.
    */
   static constexpr bool IsParameterActive( Parameter const pa ) {
      return static_cast<unsigned int>( pa ) < FieldDetails::ActiveParameters::InactiveFieldOffset;
   }

   /**
    * @brief Gives the name used for the reading of input files for the given field and index 
    * @param field_type Identifier for the material field type (Conservatives, PrimeStates or Parameters)
    * @param field_index Index of the given field 
    * @return Input name for the given field and index (empty if not specified)
    */ 
   static constexpr auto InputName( MaterialFieldType const field_type, unsigned int const field_index ) {
      switch( field_type ) {
         case MaterialFieldType::Conservatives : {
            return FieldDetails::ActiveEquations::GetInputName( field_index );
         }
         case MaterialFieldType::Parameters : {
            return FieldDetails::ActiveParameters::GetInputName( field_index );
         }
         default : { // MaterialFieldType::PrimeStates
            return FieldDetails::ActivePrimeStates::GetInputName( field_index );
         }
      }
   }

   /**
    * @brief Gives the name used in the output for the given equation.
    * @param eq Equation identifier.
    * @return Input name for the given equation (empty if not specified)
    */
   static constexpr auto InputName( Equation const eq ) {
      return FieldDetails::ActiveEquations::GetInputName( ETI( eq ) );
   }

   /**
    * @brief Gives the name used in the output for the given prime state .
    * @param ps PrimeState identifier.
    * @return Input name for the given primestate (empty if not specified)
    */
   static constexpr auto InputName( PrimeState const ps ) {
      return FieldDetails::ActivePrimeStates::GetInputName( PTI( ps ) );
   }

   /**
    * @brief Gives the name used in the output for the given parameter .
    * @param pa Parameter identifier.
    * @return Input name for the given parameter (empty if not specified)
    */
   static constexpr auto InputName( Parameter const pa ) {
      return FieldDetails::ActiveParameters::GetInputName( PTI( pa ) );
   }

   /**
    * @brief Gives the dimension/unit for the given field and index 
    * @param field_type Identifier for the material field type (Conservatives, PrimeStates or Parameters)
    * @param field_index Index of the given field 
    * @return Unit for the given field and index
    */ 
   static constexpr auto FieldUnit( MaterialFieldType const field_type, unsigned int const field_index ) {
      switch( field_type ) {
         case MaterialFieldType::Conservatives : {
            return FieldDetails::ActiveEquations::GetUnit( field_index );
         }
         case MaterialFieldType::Parameters : {
            return FieldDetails::ActiveParameters::GetUnit( field_index );
         }
         default : { // MaterialFieldType::PrimeStates
            return FieldDetails::ActivePrimeStates::GetUnit( field_index );
         }
      }
   }

   /**
    * @brief Gives the dimension/unit for the given equation.
    * @param eq Equation identifier.
    * @return Unit for the given equation 
    */
   static constexpr auto FieldUnit( Equation const eq ) {
      return FieldDetails::ActiveEquations::GetUnit( ETI( eq ) );
   }

   /**
    * @brief Gives the dimension/unit for the given prime state.
    * @param ps PrimeState identifier.
    * @return Unit for the given prime state 
    */
   static constexpr auto FieldUnit( PrimeState const ps ) {
      return FieldDetails::ActivePrimeStates::GetUnit( PTI( ps ) );
   }

   /**
    * @brief Gives the dimension/unit for the given parameter.
    * @param pa Parameter identifier.
    * @return Unit for the given parameter 
    */
   static constexpr auto FieldUnit( Parameter const pa ) {
      return FieldDetails::ActiveParameters::GetUnit( PTI( pa ) );
   }

   /**
    * @brief Gives the number of equations considered in the simulation, i.e. Euler Equations.
    * @return Number of Equations = 5 (mass, energy, X-,Y-,Z-momentum) for 3D
    * @return Number of Equations = 4 (mass, energy, X-,Y-momentum) for 2D
    * @return Number of Equations = 3 (mass, energy, X-momentum) for 1D
    */
   static constexpr unsigned int ANOE() { return FieldDetails::ActiveEquations::Count; }

   /**
    * @brief Gives the number of prime states considered in the simulation.
    * @return Number of prime states = 6 (density, pressure, temperature, X-, Y-, Z-velocity) for 3D
    * @return Number of prime states = 5 (density, pressure, temperature, X-, Y-velocity) for 2D
    * @return Number of prime states = 4 (density, pressure, temperature, X-velocity) for 1D
    */
   static constexpr unsigned int ANOP() { return FieldDetails::ActivePrimeStates::Count; }

   /**
    * @brief Gives the number of parameters considered in the simulation.
    * @return Number of active parameters
    */
   static constexpr unsigned int ANOPA() { return FieldDetails::ActiveParameters::Count; }

   /**
    * @brief Gives the number of active fields for the given field type, i.e. conservatives, prime states or parameters.
    * @param field_type The material-field type.
    * @return Number of active fields.
    */
   static constexpr unsigned int ANOF( MaterialFieldType const field_type ) {
      switch( field_type ) {
         case MaterialFieldType::Conservatives : {
            return ANOE();
         }
         case MaterialFieldType::Parameters : {
            return ANOPA();
         }
         default : { // MaterialFieldType::PrimeStates
            return ANOP();
         }
      }
   }

   /**
    * @brief Gives the set of equations which are worked on in this simulation. I. e. varies with dimension of the simulation.
    * "ASOE = Active Set of Equations".
    * @return Set of equations. E.g. Rho, Energy, X-Momentum for a 1D pure material simulation.
    */
   static constexpr auto ASOE() { return active_equations_; }

   /**
    * @brief Gives the active momentum equations, i.e. varies with dimension of the simulation. "AME = Active Momentum Equations".
    * @return Set of momentum equations.
    */
   static constexpr auto AME() { return active_momenta_; }

   /**
    * @brief Gives the set of prime states which are worked on in this simulation. I. e. varies with dimension of the simulation.
    * "ASOP = Active Set of Prime states".
    * @return Set of prime states. E.g. Rho, Pressure, X-Velocity for a 1D pure material simulation.
    */
   static constexpr auto ASOP() { return active_prime_states_; }

   /**
    * @brief Gives the active velocity prime states, i.e. varies with dimension of the simulation. "AV = Active Velocities".
    * @return Set of velocity prime states.
    */
   static constexpr auto AV() { return active_velocities_; }

   /**
    * @brief Gives the equations considered for the coarsening/refinement decision. "EWA = Equations for Wavelet-Analysis".
    * @return List of equations.
    */
   static constexpr auto EWA() { return wavelet_analysis_equations_; }

   /**
    * @brief Gives the set of parameters which are used in thi simulation "ASOPA = Active Set of PArameters".
    * @return Set of parameters.
    */
   static constexpr auto ASOPA() { return active_parameters_; }
};

using MF = MaterialFieldsDefinitions;


static_assert( std::make_pair( MF::AME()[0], MF::AV()[0] ) == std::make_pair( Equation::MomentumX, PrimeState::VelocityX ), "MF::AME()[0] and MF::AV()[0] have to consistently return MomentumX and VelocityX" );
#if DIMENSION != 1
static_assert( std::make_pair( MF::AME()[1], MF::AV()[1] ) == std::make_pair( Equation::MomentumY, PrimeState::VelocityY ), "MF::AME()[1] and MF::AV()[1] have to consistently return MomentumY and VelocityY" );
#endif
#if DIMENSION == 3
static_assert( std::make_pair( MF::AME()[2], MF::AV()[2] ) == std::make_pair( Equation::MomentumZ, PrimeState::VelocityZ ), "MF::AME()[2] and MF::AV()[2] have to consistently return MomentumZ and VelocityZ" );
#endif

#endif // FIELDS_MATERIAL_DEFINITONS_H