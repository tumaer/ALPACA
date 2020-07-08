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
#ifndef EULER_PRIME_STATE_HANDLER_H
#define EULER_PRIME_STATE_HANDLER_H

#include "prime_state_handler.h"

/**
 * @brief This class handles the conversion between conservative and prime state values for the standard Euler equations.
 * @note Also calculates temperature if PrimeState::Temperature is active.
 */
class EulerPrimeStateHandler : public PrimeStateHandler<EulerPrimeStateHandler> {
   friend PrimeStateHandler;

   template<typename PrimeStatesContainerType, typename ConservativesContainerType>
   inline void ConvertPrimeStatesToConservativesImplementation( MaterialName const& material, PrimeStatesContainerType const& prime_states_container, ConservativesContainerType& conservatives_container ) const {
      double const density                                = prime_states_container[PTI( PrimeState::Density )];
      conservatives_container[ETI( Equation::Mass )]      = density;
      conservatives_container[ETI( Equation::MomentumX )] = prime_states_container[PTI( PrimeState::VelocityX )] * density;
      if constexpr( CC::DIM() != Dimension::One ) {
         conservatives_container[ETI( Equation::MomentumY )] = prime_states_container[PTI( PrimeState::VelocityY )] * density;
      }
      if constexpr( CC::DIM() == Dimension::Three ) {
         conservatives_container[ETI( Equation::MomentumZ )] = prime_states_container[PTI( PrimeState::VelocityZ )] * density;
      }
      conservatives_container[ETI( Equation::Energy )] = material_manager_.GetMaterial( material ).GetEquationOfState().GetEnergy( density, conservatives_container[ETI( Equation::MomentumX )], CC::DIM() != Dimension::One ? conservatives_container[ETI( Equation::MomentumY )] : 0.0, CC::DIM() == Dimension::Three ? conservatives_container[ETI( Equation::MomentumZ )] : 0.0, prime_states_container[PTI( PrimeState::Pressure )] );
   }

   template<typename ConservativesContainerType, typename PrimeStatesContainerType>
   inline void ConvertConservativesToPrimeStatesImplementation( MaterialName const& material, ConservativesContainerType const& conservatives_container, PrimeStatesContainerType& prime_states_container ) const {
      double const one_density                             = 1.0 / conservatives_container[ETI( Equation::Mass )];
      prime_states_container[PTI( PrimeState::Density )]   = conservatives_container[ETI( Equation::Mass )];
      prime_states_container[PTI( PrimeState::VelocityX )] = conservatives_container[ETI( Equation::MomentumX )] * one_density;
      if constexpr( CC::DIM() != Dimension::One ) {
         prime_states_container[PTI( PrimeState::VelocityY )] = conservatives_container[ETI( Equation::MomentumY )] * one_density;
      }
      if constexpr( CC::DIM() == Dimension::Three ) {
         prime_states_container[PTI( PrimeState::VelocityZ )] = conservatives_container[ETI( Equation::MomentumZ )] * one_density;
      }
      prime_states_container[PTI( PrimeState::Pressure )] = material_manager_.GetMaterial( material ).GetEquationOfState().GetPressure( conservatives_container[ETI( Equation::Mass )], conservatives_container[ETI( Equation::MomentumX )], CC::DIM() != Dimension::One ? conservatives_container[ETI( Equation::MomentumY )] : 0.0, CC::DIM() == Dimension::Three ? conservatives_container[ETI( Equation::MomentumZ )] : 0.0, conservatives_container[ETI( Equation::Energy )] );
      if constexpr( MF::IsPrimeStateActive( PrimeState::Temperature ) ) {
         // only calculate temperature if the prime state is activated
         prime_states_container[PTI( PrimeState::Temperature )] = material_manager_.GetMaterial( material ).GetEquationOfState().GetTemperature( conservatives_container[ETI( Equation::Mass )], conservatives_container[ETI( Equation::MomentumX )], CC::DIM() != Dimension::One ? conservatives_container[ETI( Equation::MomentumY )] : 0.0, CC::DIM() == Dimension::Three ? conservatives_container[ETI( Equation::MomentumZ )] : 0.0, conservatives_container[ETI( Equation::Energy )] );
      }
   }

public:
   ~EulerPrimeStateHandler() = default;
   /**
    * @brief Construct a new EulerPrimeStateHandler using an existing MaterialManager instance.
    * @param material_manager The MaterialManager object to be used for conversion.
    */
   explicit EulerPrimeStateHandler( MaterialManager const& material_manager ) : PrimeStateHandler( material_manager ) {}
   EulerPrimeStateHandler( EulerPrimeStateHandler const& ) = delete;
   EulerPrimeStateHandler& operator=( EulerPrimeStateHandler const& ) = delete;
   EulerPrimeStateHandler( EulerPrimeStateHandler&& )                 = delete;
   EulerPrimeStateHandler& operator=( EulerPrimeStateHandler&& ) = delete;
};

#endif// EULER_PRIME_STATE_HANDLER_H
