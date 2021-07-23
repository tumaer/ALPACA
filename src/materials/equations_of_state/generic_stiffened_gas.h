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
#ifndef GENERIC_STIFFENED_GAS_H
#define GENERIC_STIFFENED_GAS_H

#include "utilities/mathematical_functions.h"
#include <cmath>

namespace GenericStiffenedGas {
   /**
    * @brief Computes pressure from inputs as -gamma*pi + (gamma - 1) * (E  - 0.5 * rho * ||v^2||)
    * @tparam safe decision flag whether to take care of devision by 0
    * @param density .
    * @param momentum_x .
    * @param momentum_y .
    * @param momentum_z .
    * @param energy .
    * @return Pressure according to stiffened-gas equation of state.
    */
   template<bool safe>
   constexpr double CalculatePressure( const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy, const double gamma, const double pi ) {
      double const q = 0.5 * DimensionAwareConsistencyManagedSum( momentum_x * momentum_x, momentum_y * momentum_y, momentum_z * momentum_z );
      if constexpr( safe ) {
         double const pressure = -gamma * pi + ( gamma - 1.0 ) * ( energy - q / std::max( density, std::numeric_limits<double>::epsilon() ) );
         return std::max( pressure, -pi + std::numeric_limits<double>::epsilon() );
      }
      return -gamma * pi + ( gamma - 1.0 ) * ( energy - q / density );
   }

   /**
   * @brief Computes energy from inputs as (p + gamma * pi) / (gamma - 1) + 0.5 * rho * ||v^2||
   * @param density .
   * @param velocity_x .
   * @param velocity_y .
   * @param velocity_z .
   * @param pressure .
   * @return Energy according to stiffened-gas equation of state.
   */
   constexpr double CalculateEnergy( const double density, const double velocity_x, const double velocity_y, const double velocity_z, const double pressure, double const gamma, double const pi ) {
      return ( pressure + gamma * pi ) / ( gamma - 1.0 ) + ( 0.5 * DimensionAwareConsistencyManagedSum( velocity_x * velocity_x, velocity_y * velocity_y, velocity_z * velocity_z ) * density );
   }

   /**
   * @brief Computes speed of sound from inputs as sqrt(gamma * (p + pi)) / rho
   * @tparam safe decision flag whether to take care of devision by 0 and negative root
   * @param density .
   * @param pressure .
   * @return Speed of sound according to stiffened-gas equation of state.
   */
   template<bool safe>
   constexpr double CalculateSpeedOfSound( const double density, const double pressure, const double gamma, const double pi ) {
      if constexpr( safe ) {
         const double speed_of_sound_squared = gamma * ( pressure + pi ) / std::max( density, std::numeric_limits<double>::epsilon() );
         return std::sqrt( std::max( speed_of_sound_squared, std::numeric_limits<double>::epsilon() ) );
      } else {
         return std::sqrt( gamma * ( pressure + pi ) / density );
      }
   }

}// namespace GenericStiffenedGas

#endif//GENERIC_STIFFENED_GAS_H
