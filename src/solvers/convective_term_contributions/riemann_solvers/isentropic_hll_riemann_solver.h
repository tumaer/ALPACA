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

#ifndef ISENTROPIC_HLL_RIEMANN_SOLVER_H
#define ISENTROPIC_HLL_RIEMANN_SOLVER_H

#include "solvers/convective_term_contributions/riemann_solvers/riemann_solver.h"
#include "block_definitions/block.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief Discretization of the Riemann solver using the HLL procedure derived from \cite Toro2009, chapter 10.3 for isentropic cases.
 */
class IsentropicHllRiemannSolver : public RiemannSolver<IsentropicHllRiemannSolver> {

   friend RiemannSolver;

   /**
    * @brief Computes an approximate solution of a Riemann problem using the HLL procedure.
    * @param material Material information of the phase under consideration.
    * @param conservatives_left The initial left conservative states.
    * @param conservatives_right The initial right conservative states.
    * @param prime_state_left The initial left primitive states.
    * @param prime_state_right The initial right primitive states.
    * @tparam DIR Indicates which spatial direction is to be computed.
    * @note Hotpath function.
    */
   template<Direction DIR>
   std::array<double, MF::ANOE()> SolveRiemannProblemImplementation( MaterialName const material,
                                                                     std::array<double, MF::ANOE()> const& conservatives_left,
                                                                     std::array<double, MF::ANOE()> const& conservatives_right,
                                                                     std::array<double, MF::ANOP()> const& prime_state_left,
                                                                     std::array<double, MF::ANOP()> const& prime_state_right ) const {
      std::array<double, MF::ANOE()> flux_left;
      std::array<double, MF::ANOE()> flux_right;
      std::array<double, MF::ANOE()> fluxes;

      constexpr unsigned int principal_momentum_index = ETI( MF::AME()[DTI( DIR )] );
      constexpr unsigned int principal_velocity_index = PTI( MF::AV()[DTI( DIR )] );

      // For Toro signal speeds
      double const gamma = material_manager_.GetMaterial( material ).GetEquationOfState().Gamma();

      // Compute pressure, velocity and speed of sound for both cells for reconstructed values
      double const one_density_left     = 1.0 / conservatives_left[ETI( Equation::Mass )];
      double const one_density_right    = 1.0 / conservatives_right[ETI( Equation::Mass )];
      double const pressure_left        = prime_state_left[PTI( PrimeState::Pressure )];
      double const pressure_right       = prime_state_right[PTI( PrimeState::Pressure )];
      double const velocity_left        = prime_state_left[principal_velocity_index];
      double const velocity_right       = prime_state_right[principal_velocity_index];
      double const speed_of_sound_left  = material_manager_.GetMaterial( material ).GetEquationOfState().SpeedOfSound( conservatives_left[ETI( Equation::Mass )], pressure_left );
      double const speed_of_sound_right = material_manager_.GetMaterial( material ).GetEquationOfState().SpeedOfSound( conservatives_right[ETI( Equation::Mass )], pressure_right );

      // Calculation of signal speeds
      auto const [wave_speed_left_simple, wave_speed_right_simple] = CalculateSignalSpeed( conservatives_left[ETI( Equation::Mass )], conservatives_right[ETI( Equation::Mass )],
                                                                                           velocity_left, velocity_right,
                                                                                           pressure_left, pressure_right,
                                                                                           speed_of_sound_left, speed_of_sound_right,
                                                                                           gamma );

      double const wave_speed_left  = std::min( wave_speed_left_simple, 0.0 );
      double const wave_speed_right = std::max( wave_speed_right_simple, 0.0 );

      // Calculation of left and right flux
      flux_left[ETI( Equation::Mass )]     = conservatives_left[principal_momentum_index];
      flux_left[principal_momentum_index]  = ( ( conservatives_left[principal_momentum_index] * conservatives_left[principal_momentum_index] ) * one_density_left ) + pressure_left;
      flux_right[ETI( Equation::Mass )]    = conservatives_right[principal_momentum_index];
      flux_right[principal_momentum_index] = ( ( conservatives_right[principal_momentum_index] * conservatives_right[principal_momentum_index] ) * one_density_right ) + pressure_right;

      // minor momenta
      for( unsigned int d = 0; d < DTI( CC::DIM() ) - 1; ++d ) {
         // get the index of this minor momentum
         unsigned int const minor_momentum_index = ETI( MF::AME()[DTI( GetMinorDirection<DIR>( d ) )] );
         flux_left[minor_momentum_index]         = velocity_left * conservatives_left[minor_momentum_index];
         flux_right[minor_momentum_index]        = velocity_right * conservatives_right[minor_momentum_index];
      }

      for( unsigned int n = 0; n < MF::ANOE(); ++n ) {
         fluxes[n] = ( ( ( wave_speed_right * flux_left[n] ) - ( wave_speed_left * flux_right[n] ) ) +
                       ( ( wave_speed_right * wave_speed_left ) * ( conservatives_right[n] - conservatives_left[n] ) ) ) /
                     ( wave_speed_right - wave_speed_left );
      }
      return fluxes;
   }

public:
   IsentropicHllRiemannSolver() = delete;
   explicit IsentropicHllRiemannSolver( MaterialManager const& material_manager, EigenDecomposition const& eigendecomposition_calculator ) : RiemannSolver( material_manager, eigendecomposition_calculator ) {}
   ~IsentropicHllRiemannSolver()                                   = default;
   IsentropicHllRiemannSolver( IsentropicHllRiemannSolver const& ) = delete;
   IsentropicHllRiemannSolver& operator=( IsentropicHllRiemannSolver const& ) = delete;
   IsentropicHllRiemannSolver( IsentropicHllRiemannSolver&& )                 = delete;
   IsentropicHllRiemannSolver& operator=( IsentropicHllRiemannSolver&& ) = delete;
};

#endif// HLL_RIEMANN_SOLVER_H
