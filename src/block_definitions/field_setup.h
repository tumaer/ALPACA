
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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#ifndef FIELD_SETUP_H
#define FIELD_SETUP_H

#include "utilities/helper_functions.h"
#include "user_specifications/equation_settings.h"
#include "enums/unit_type.h"
#include "block_definitions/field_enums.h"
#include <string_view>

namespace {
   /**
 * @brief Dummy to allow empty Arrays - as given for demonstartion purposes only (!) - in custom field definitions.
 * Will not be used if custom equations are provided.
 */
   constexpr std::array<void*, 0> MakeArray() { return {}; }
}// namespace

namespace FieldSetup {

   template<typename FieldEnum>
   struct FieldInfo {
      FieldEnum const Field;
      UnitType const Unit;
      std::string_view const Name = "";
      // TODO Rename members and constructo inputs. Constructo given, because it allows autoamtic template deduction.
      constexpr FieldInfo( FieldEnum const fe, UnitType const ut, std::string_view const in ) : Field( fe ), Unit( ut ), Name( in ) {}
   };

   namespace Custom {
      constexpr auto equations = MakeArray(
            // clang-format off
            //                         conservative                  unit             name
            //   FieldInfo{ EquationPool::Example,           UnitType::Example,      "example" }
            // , FieldInfo{ EquationPool::Example2,          UnitType::Example2,     "example2" }
            // clang-format on
      );
      constexpr auto primes = MakeArray(
            // clang-format off
            //              prime state                     unit                       name
            //  FieldInfo{ PrimeStatePool::Example,         UnitType::Example,         "example" }
            //, FieldInfo{ PrimeStatePool::Example2,        UnitType::Example2,        "" }
            // clang-format on
      );
   }// namespace Custom

   namespace Isentropic {
      constexpr auto equations = MakeArray(
            // clang-format off
            //           conservative           unit                input name
            FieldInfo{ EquationPool::Mass,      UnitType::Density,  "" },
            FieldInfo{ EquationPool::MomentumX, UnitType::Momentum, "" }
#if DIMENSION != 1
           ,FieldInfo{ EquationPool::MomentumY, UnitType::Momentum, "" }
#endif
#if DIMENSION == 3
           ,FieldInfo{ EquationPool::MomentumZ, UnitType::Momentum, "" }
#endif
            // clang-format on
      );

      constexpr auto primes = MakeArray(
            // clang-format off
            //           prime state              unit                name
            FieldInfo{ PrimeStatePool::Density,   UnitType::Density,  "density" },
            FieldInfo{ PrimeStatePool::Pressure,  UnitType::Pressure, "" },
            FieldInfo{ PrimeStatePool::VelocityX, UnitType::Velocity, "velocityX" }
#if DIMENSION != 1
           ,FieldInfo{ PrimeStatePool::VelocityY, UnitType::Velocity, "velocityY" }
#endif
#if DIMENSION == 3
           ,FieldInfo{ PrimeStatePool::VelocityZ, UnitType::Velocity, "velocityZ" }
#endif
            // clang-format on
      );
   }// namespace Isentropic

   namespace Euler {
      constexpr auto equations = MakeArray(
            // clang-format off
            //           conservative           unit                name
            FieldInfo{ EquationPool::Mass,      UnitType::Density,  "" },
            FieldInfo{ EquationPool::Energy,    UnitType::Energy,   "" },
            FieldInfo{ EquationPool::MomentumX, UnitType::Momentum, "" }
#if DIMENSION != 1
           ,FieldInfo{ EquationPool::MomentumY, UnitType::Momentum, "" }
#endif
#if DIMENSION == 3
           ,FieldInfo{ EquationPool::MomentumZ, UnitType::Momentum, "" }
#endif
            // clang-format on
      );

      constexpr auto primes = MakeArray(
            // clang-format off
            //         prime state                unit                  name
            FieldInfo{ PrimeStatePool::Density,   UnitType::Density,  "density" },
            FieldInfo{ PrimeStatePool::Pressure,  UnitType::Pressure, "pressure" },
            FieldInfo{ PrimeStatePool::VelocityX, UnitType::Velocity, "velocityX" }
#if DIMENSION != 1
           ,FieldInfo{ PrimeStatePool::VelocityY, UnitType::Velocity, "velocityY" }
#endif
#if DIMENSION == 3
           ,FieldInfo{ PrimeStatePool::VelocityZ, UnitType::Velocity, "velocityZ" }
#endif
            // clang-format on
      );
   }// namespace Euler

   namespace NavierStokes {
      constexpr auto equations = MakeArray(
            // clang-format off
            //         conservative             unit                name
            FieldInfo{ EquationPool::Mass,      UnitType::Density,  "" },
            FieldInfo{ EquationPool::Energy,    UnitType::Energy,   "" },
            FieldInfo{ EquationPool::MomentumX, UnitType::Momentum, "" }
#if DIMENSION != 1
           ,FieldInfo{ EquationPool::MomentumY, UnitType::Momentum, "" }
#endif
#if DIMENSION == 3
           ,FieldInfo{ EquationPool::MomentumZ, UnitType::Momentum, "" }
#endif
            // clang-format on
      );

      constexpr auto primes = MakeArray(
            // clang-format off
            //         prime state                      unit               name
            FieldInfo{ PrimeStatePool::Density,     UnitType::Density,     "density" },
            FieldInfo{ PrimeStatePool::Pressure,    UnitType::Pressure,    "pressure" },
            FieldInfo{ PrimeStatePool::VelocityX,   UnitType::Velocity,    "velocityX" }
#if DIMENSION != 1
           ,FieldInfo{ PrimeStatePool::VelocityY,   UnitType::Velocity,    "velocityY" }
#endif
#if DIMENSION == 3
           ,FieldInfo{ PrimeStatePool::VelocityZ,   UnitType::Velocity,    "velocityZ" }
#endif
           ,FieldInfo{ PrimeStatePool::Temperature, UnitType::Temperature, "" } );
      // clang-format on
   }// namespace NavierStokes

   template<EquationSet>
   constexpr auto EquationDefinition();

   template<EquationSet>
   constexpr auto PrimeStateDefinition();

   template<>
   constexpr auto EquationDefinition<EquationSet::Isentropic>() {
      return Isentropic::equations;
   }

   template<>
   constexpr auto PrimeStateDefinition<EquationSet::Isentropic>() {
      return Isentropic::primes;
   }

   template<>
   constexpr auto EquationDefinition<EquationSet::Euler>() {
      return Euler::equations;
   }

   template<>
   constexpr auto PrimeStateDefinition<EquationSet::Euler>() {
      return Euler::primes;
   }

   template<>
   constexpr auto EquationDefinition<EquationSet::NavierStokes>() {
      return NavierStokes::equations;
   }

   template<>
   constexpr auto PrimeStateDefinition<EquationSet::NavierStokes>() {
      return NavierStokes::primes;
   }

   template<>
   constexpr auto EquationDefinition<EquationSet::Custom>() {
      return Custom::equations;
   }

   template<>
   constexpr auto PrimeStateDefinition<EquationSet::Custom>() {
      return Custom::primes;
   }

   template<EquationSet>
   constexpr auto WaveletEquations() {
      std::array<EquationPool, 2> equations = { EquationPool::Mass, EquationPool::Energy };
      return equations;
   }

   template<>
   constexpr auto WaveletEquations<EquationSet::Isentropic>() {
      std::array<EquationPool, 1> equations = { EquationPool::Mass };
      return equations;
   }

   namespace Interface {
      namespace State {
         constexpr auto custom_quantities = MakeArray(
               // clang-format off
               //           interface quantity                          unit                  name
               //   FieldInfo{ InterfaceStatePool::Example,             UnitType::Example,    "example" }
               // , FieldInfo{ InterfaceStatePool::Example2,            UnitType::Example2,   "example2" }
               // clang-format on
         );
         constexpr auto default_quantities = MakeArray(
               // clang-format off
               //           interface quantity                    unit                  name
               FieldInfo{ InterfaceStatePool::Velocity,         UnitType::Velocity, "" },
               FieldInfo{ InterfaceStatePool::PressurePositive, UnitType::Pressure, "" },
               FieldInfo{ InterfaceStatePool::PressureNegative, UnitType::Pressure, "" } );
         // clang-format on
      }// namespace State
      namespace Description {
         static constexpr auto custom_quantities = MakeArray(
               // clang-format off
               //           interface description                               unit                  name
               // FieldInfo{ InterfaceDescriptionPool::Example,             UnitType::Example,      "example" }
               // clang-format on
         );

         static constexpr auto default_quantities = MakeArray(
               // clang-format off
               //           interface description                    unit                name
               FieldInfo{ InterfaceDescriptionPool::Levelset,       UnitType::Unitless, "phi" },
               FieldInfo{ InterfaceDescriptionPool::VolumeFraction, UnitType::Unitless, "" } );
         // clang-format on
      }// namespace Description

      namespace Parameter {
         static constexpr auto custom_quantities = MakeArray(
               // clang-format off
               //           interface parameter                   unit                 name
               //FieldInfo{ InterfaceParameterPool::Example,      UnitType::Example,   "example" }
               // clang-format on
         );

         static constexpr auto default_quantities = MakeArray(
               // clang-format off
               //         interface parameter                                unit                                 name
               FieldInfo{ InterfaceParameterPool::SurfaceTensionCoefficient, UnitType::SurfaceTensionCoefficient, "" } );
         // clang-format on
      }// namespace Parameter

   }// namespace Interface

   template<InterfaceSet>
   constexpr auto InterFaceStateDefinition();

   template<>
   constexpr auto InterFaceStateDefinition<InterfaceSet::Default>() {
      return Interface::State::default_quantities;
   }

   template<>
   constexpr auto InterFaceStateDefinition<InterfaceSet::Custom>() {
      return Interface::State::custom_quantities;
   }

   template<InterfaceSet>
   constexpr auto InterFaceDescriptionDefinition();

   template<>
   constexpr auto InterFaceDescriptionDefinition<InterfaceSet::Default>() {
      return Interface::Description::default_quantities;
   }

   template<>
   constexpr auto InterFaceDescriptionDefinition<InterfaceSet::Custom>() {
      return Interface::Description::custom_quantities;
   }

   template<InterfaceSet>
   constexpr auto InterFaceParameterDefinition();

   template<>
   constexpr auto InterFaceParameterDefinition<InterfaceSet::Default>() {
      return Interface::Parameter::default_quantities;
   }

   template<>
   constexpr auto InterFaceParameterDefinition<InterfaceSet::Custom>() {
      return Interface::Parameter::custom_quantities;
   }

   namespace Parameters {
      static constexpr auto custom_quantities = MakeArray(
            // clang-format off
            //         parameter                                      unit                  input_name
            //FieldInfo{ ParameterPool::Example,                    UnitType::Example,      "example" }
            // clang-format on
      );

      static constexpr auto default_quantities = MakeArray(
            // clang-format off
            //         parameter                                      unit                  input_name
            FieldInfo{ ParameterPool::ShearViscosity, UnitType::Viscosity, "" },
            FieldInfo{ ParameterPool::ThermalConductivity, UnitType::ThermalConductivity, "" } );
      // clang-format on
   }// namespace Parameters

   template<ParameterSet>
   constexpr auto ParameterDefinition();

   template<>
   constexpr auto ParameterDefinition<ParameterSet::Default>() {
      return Parameters::default_quantities;
   }

   template<>
   constexpr auto ParameterDefinition<ParameterSet::Custom>() {
      return Parameters::custom_quantities;
   }
}// namespace FieldSetup

#endif//FIELD_SETUP_H
