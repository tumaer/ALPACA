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
#include "levelset/geometry/geometry_calculator.h"
#include "capillary_pressure_calculator.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "utilities/mathematical_functions.h"
#include "stencils/stencil_utilities.h"
#include "enums/interface_tag_definition.h"

/**
 * @brief      The default constructor of the class.
 *
 * @param[in]  material_manager  The material manager containing information about the materials used in the simulation. Necessary to get the
 * surface tension coefficient.
 */
CapillaryPressureCalculator::CapillaryPressureCalculator(MaterialManager const& material_manager) :
   surface_tension_coefficient_( CC::CapillaryForcesActive() ? 
                                 material_manager.GetMaterialPairing( MaterialSignCapsule::NegativeMaterial(), 
                                                                      MaterialSignCapsule::PositiveMaterial()).GetSurfaceTensionCoefficient() : 0.0 ) {
   // Empty besides initializer list.
}

/**
 * @brief Calculates the capillary pressure in cut cells. Stores the result in a buffer.
 *
 * @param node                 The node for which the capillary pressure is calculated.
 * @param pressure_difference  The buffer in which the capillary pressure is saved. Indirect return parameter.
 */
void CapillaryPressureCalculator::ComputePressureDifference( Node const& node, double (&pressure_difference)[CC::TCX()][CC::TCY()][CC::TCZ()] ) const {

   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   /**
    * Use the pressure_difference buffer in the next function call to store the curvature. Then, in the following step
    * multiplication with the surface_tension_coefficient_ gives the pressure difference due to capillary forces.
    */
   ComputeInterfaceCurvature( node, pressure_difference );

   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell ) ) {
               pressure_difference[i][j][k] *= surface_tension_coefficient_;
            }
         }
      }
   }

}
