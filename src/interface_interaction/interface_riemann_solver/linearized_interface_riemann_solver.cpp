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
#include "linearized_interface_riemann_solver.h"
#include <limits>

/**
 * @brief Constructor for the LinearizedInterfaceRiemannSolver.
 * @param material_manager See base class.
 */
LinearizedInterfaceRiemannSolver::LinearizedInterfaceRiemannSolver( MaterialManager const& material_manager ) : InterfaceRiemannSolver( material_manager ) {
   // Empty besides call of base class constructor.
}

/**
 * @brief Solves the Riemann problem at the interface linearized. See \cite Saurel2003 for details.
 * @param rho_left Density of the left fluid.
 * @param p_left Pressure of the left fluid.
 * @param velocity_normal_left Velocity normal to the interface of the left fluid.
 * @param material_left Material of the left fluid.
 * @param rho_right Density of the right fluid.
 * @param p_right Pressure of the right fluid.
 * @param velocity_normal_right Velocity normal to the interface of the right fluid.
 * @param material_right Material of the right fluid.
 * @param delta_p Pressure jump due to capillarity.
 * @return An array that contains following information in the given order: interface_velocity, interface_pressure_positive, interface_pressure_negative.
 */
std::array<double, 3> LinearizedInterfaceRiemannSolver::SolveInterfaceRiemannProblemImplementation( double const rho_left, double const p_left, double const velocity_normal_left, MaterialName const material_left,
                                                                                                    double const rho_right, double const p_right, double const velocity_normal_right, MaterialName const material_right,
                                                                                                    double const delta_p ) const {
   double const c_left         = material_manager_.GetMaterial( material_left ).GetEquationOfState().SpeedOfSound( rho_left, p_left );
   double const impedance_left = rho_left * c_left;

   double const c_right         = material_manager_.GetMaterial( material_right ).GetEquationOfState().SpeedOfSound( rho_right, p_right );
   double const impedance_right = rho_right * c_right;

   double const inverse_impedance_sum = 1.0 / std::max( ( impedance_left + impedance_right ), std::numeric_limits<double>::epsilon() );

   double const interface_velocity = ( impedance_left * velocity_normal_left +
                                       impedance_right * velocity_normal_right +
                                       p_left - p_right - delta_p ) *
                                     inverse_impedance_sum;

   double const interface_pressure_positive = ( impedance_left * p_right +
                                                impedance_right * ( p_left - delta_p ) +
                                                impedance_left * impedance_right * ( velocity_normal_left - velocity_normal_right ) ) *
                                              inverse_impedance_sum;

   double const interface_pressure_negative = ( impedance_left * ( p_right + delta_p ) +
                                                impedance_right * p_left +
                                                impedance_left * impedance_right * ( velocity_normal_left - velocity_normal_right ) ) *
                                              inverse_impedance_sum;

   return { interface_velocity, interface_pressure_positive, interface_pressure_negative };
}
