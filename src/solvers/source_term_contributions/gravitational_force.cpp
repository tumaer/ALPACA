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
#include "gravitational_force.h"

#include "utilities/index_transformations.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief The default constructor of the class.
 * @param gravity The vector containing the components of the gravitational force.
 */
GravitationalForce::GravitationalForce( std::array<double, 3> const gravity ) :
   gravity_(gravity)
{
   // Empty besides initializer list
}

/**
 * @brief Computes increments for cell averages due to gravity.
 * @param block Block of the considered phase.
 * @param gravity_forces Reference to array of volume forces increments to be filled here (indirect return parameter).
 */
void GravitationalForce::ComputeForces( Block const& block, double (&gravity_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()] ) const {

   Conservatives const& conservatives = block.GetAverageBuffer();
   double const(&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer( Equation::Mass );

   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {

            std::array<unsigned int, 3> const indices = { BIT::I2TX(i), BIT::I2TY(j), BIT::I2TZ(k) };

            // Add up to volume forces
            gravity_forces[ETI(Equation::Mass)][i][j][k] += 0.0;
            gravity_forces[ETI(Equation::Energy) ][i][j][k] += DimensionAwareConsistencyManagedSum( gravity_[0] * conservatives[Equation::MomentumX][indices[0]][indices[1]][indices[2]]
                                                       , CC::DIM() != Dimension::One   ? gravity_[1] * conservatives[Equation::MomentumY][indices[0]][indices[1]][indices[2]] : 0.0
                                                       , CC::DIM() == Dimension::Three ? gravity_[2] * conservatives[Equation::MomentumZ][indices[0]][indices[1]][indices[2]] : 0.0 );

            gravity_forces[ETI(Equation::MomentumX)][i][j][k] += gravity_[0] * density[indices[0]][indices[1]][indices[2]];
            if( CC::DIM() != Dimension::One ) {
               gravity_forces[ETI(Equation::MomentumY)][i][j][k] += gravity_[1] * density[indices[0]][indices[1]][indices[2]];
            }
            if( CC::DIM() == Dimension::Three ) {
               gravity_forces[ETI(Equation::MomentumZ)][i][j][k] += gravity_[2] * density[indices[0]][indices[1]][indices[2]];
            }
         }
      }
   }
}
