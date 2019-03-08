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
#include "interface_state_calculator.h"
#include "levelset/multi_phase_manager/two_phase_manager.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief A constructor for the InterfaceStateCalculator.
 * @param material_manager A material manager giving information about the materials present in the simulation.
 */
InterfaceStateCalculator::InterfaceStateCalculator(MaterialManager const& material_manager) :
   interface_velocity_pressure_calculator_(material_manager),
   capillary_pressure_calculator_(material_manager) {
  // Empty besides initializer list.
}

/**
 * @brief Computes the interface quantities for a node.
 * @param node The node for which the interface quantities are calculated.
 * @note This function should only be called for multi-phase nodes
 * that have a level-set block. Sanity check whether the input is a multi-phase node is not done!
 */
void InterfaceStateCalculator::ObtainInterfaceQuantities(Node& node) const {

  double pressure_difference[CC::TCX()][CC::TCY()][CC::TCZ()];
  for( unsigned int i = 0; i < CC::TCX(); ++i ){
    for( unsigned int j = 0; j < CC::TCY(); ++j ){
      for( unsigned int k = 0; k < CC::TCZ(); ++k ){
        pressure_difference[i][j][k] = 0.0;
      }
    }
  }

  if( CC::CapillaryForcesActive() && CC::DIM() != Dimension::One ) {
    capillary_pressure_calculator_.ComputePressureDifference( node, pressure_difference );
  }
  interface_velocity_pressure_calculator_.FillInterfaceVelocityAndPressureBuffer( node, pressure_difference );

}
