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
#ifndef RUNGE_KUTTA_2_TVD_H
#define RUNGE_KUTTA_2_TVD_H

#include "time_integrator.h"

/**
 * @brief The RungeKutta2TVD class integrates in time using a total variation diminishing two-step Runge-Kutta method. On paper the equations are
 *        I)    u^*     =       u^n +           +       dt * f(u^n)
 *        II)   u^(n+1) = 0.5 * u^n + 0.5 * u^* + 0.5 * dt * f(u^*)
 *        To reduce the memory footprint the equations are reordered to allow the integration using two buffers only.
 */
class RungeKutta2TVD : public TimeIntegrator<RungeKutta2TVD> {

   friend TimeIntegrator;

   static constexpr unsigned int number_of_stages_ = 2;

   static constexpr std::array<double, number_of_stages_> timestep_multiplier_jump_conservatives_ = { 0.5, // first stage
                                                                                                      0.5  // second stage
                                                                                                    };
                                                                           
   static constexpr std::array<double, number_of_stages_> timestep_multiplier_conservatives_      = { 1.0, // first stage
                                                                                                      0.5  // second stage 
                                                                                                    };

   static constexpr std::array<std::array<double, 2>, number_of_stages_ - 1> buffer_multiplier_       = {{
                                                                                                       {0.5, 0.5 }  // second stage
                                                                                                      }};

public:
   RungeKutta2TVD() = delete;
   ~RungeKutta2TVD() = default;
   RungeKutta2TVD( RungeKutta2TVD const& ) = delete;
   RungeKutta2TVD& operator=( RungeKutta2TVD const& ) = delete;
   RungeKutta2TVD( RungeKutta2TVD&& ) = delete;
   RungeKutta2TVD& operator=( RungeKutta2TVD&& ) = delete;

   /**
     * @brief Constructor.
     * @param start_time Time when the simulation should start.
     */
   explicit RungeKutta2TVD( double const start_time = 0.0 ) : TimeIntegrator(start_time) {}

};

#endif // RUNGE_KUTTA_2_TVD_H
