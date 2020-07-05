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
#include "interface_velocity_pressure_calculator.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "levelset/geometry/geometry_calculator.h"
#include "enums/interface_tag_definition.h"

/**
 * @brief      Default constructor for the InterfaceVelocityPressureCalculator class. Initializes the members.
 *
 * @param[in]  material_manager  The material manager containing information about the materials involved in the simulation.
 */
InterfaceVelocityPressureCalculator::InterfaceVelocityPressureCalculator( MaterialManager const& material_manager ) :
   interface_riemann_solver_( material_manager ) {
   // Empty besides initializer list
}

/**
 * @brief      Fills the interface velocity and pressure buffer of the levelset block of a given node. This function may consider a pressure jump due to capillary forces
 * at the interface.
 *
 * @param      node                 The node for which the interface velocity and pressure are calculated.
 * @param      pressure_difference  The pressure difference induced by capillary forces.
 */
void InterfaceVelocityPressureCalculator::FillInterfaceVelocityAndPressureBuffer( Node& node, double (&pressure_difference)[CC::TCX()][CC::TCY()][CC::TCZ()] ) const {

   double const (&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer(InterfaceDescription::Levelset);
   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();

   double (&interface_velocity)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetInterfaceStateBuffer(InterfaceState::Velocity);
   double (&interface_pressure_positive)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetInterfaceStateBuffer(InterfaceState::PressurePositive);
   double (&interface_pressure_negative)[CC::TCX()][CC::TCY()][CC::TCZ()] = CC::CapillaryForcesActive() ?
                                                                            node.GetInterfaceBlock().GetInterfaceStateBuffer(InterfaceState::PressureNegative) :
                                                                            node.GetInterfaceBlock().GetInterfaceStateBuffer(InterfaceState::PressurePositive);

   // for describing positive material as right, and negative material als left, see cited paper of Luo
   MaterialName const material_left  = MaterialSignCapsule::NegativeMaterial();
   MaterialName const material_right = MaterialSignCapsule::PositiveMaterial();

   PrimeStates const& left_prime_states = node.GetPhaseByMaterial( material_left ).GetPrimeStateBuffer();
   PrimeStates const& right_prime_states = node.GetPhaseByMaterial( material_right ).GetPrimeStateBuffer();

   double velocity_normal_left = 0.0;
   double velocity_normal_right = 0.0;

   for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
         for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            if(std::abs(interface_tags[i][j][k]) <= ITTI(IT::ExtensionBand)){

               std::array<double, 3> const normal = GetNormal(levelset_reinitialized, i, j, k);

               velocity_normal_left = left_prime_states[PrimeState::VelocityX][i][j][k] * normal[0];
               velocity_normal_left += CC::DIM() != Dimension::One   ? left_prime_states[PrimeState::VelocityY][i][j][k] * normal[1] : 0.0;
               velocity_normal_left += CC::DIM() == Dimension::Three ? left_prime_states[PrimeState::VelocityZ][i][j][k] * normal[2] : 0.0;

               velocity_normal_right = right_prime_states[PrimeState::VelocityX][i][j][k] * normal[0];
               velocity_normal_right += CC::DIM() != Dimension::One   ? right_prime_states[PrimeState::VelocityY][i][j][k] * normal[1] : 0.0;
               velocity_normal_right += CC::DIM() == Dimension::Three ? right_prime_states[PrimeState::VelocityZ][i][j][k] * normal[2] : 0.0;

               std::array<double, 3> const interface_states = interface_riemann_solver_.SolveInterfaceRiemannProblem( left_prime_states[PrimeState::Density][i][j][k], left_prime_states[PrimeState::Pressure][i][j][k], velocity_normal_left, material_left,
                  right_prime_states[PrimeState::Density][i][j][k], right_prime_states[PrimeState::Pressure][i][j][k], velocity_normal_right, material_right,
                  pressure_difference[i][j][k] );

               interface_velocity[i][j][k] = interface_states[0];
               interface_pressure_positive[i][j][k] = interface_states[1];
               if constexpr( CC::CapillaryForcesActive() ) {
                  interface_pressure_negative[i][j][k] = interface_states[2];
               }
            } //if
         } //k
      } //j
   } //i
}
