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
#include "axisymmetric_fluxes.h"

/**
 * @brief Computes source terms for axisymmetric simulations. Terms are set accoring to \cite Adami2016 .
 * @param block Block of the considered phase.
 * @param volume_forces Reference to array of volume forces increments to be filled here (indirect return parameter).
 */
void AxisymmetricFluxes::ComputeAxisymmetricContributions( Block const& block, double (&volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()],
   double const cell_size, double const x_block_coordinate ) const {

   double const (&velocity_x)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetPrimeStateBuffer(PrimeState::VelocityX);
   double const   (&pressure)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetPrimeStateBuffer(PrimeState::Pressure);
   double const (&momentum_x)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(Equation::MomentumX);
   // direct use of y-momentum buffer allowed since axisymmetric is only used in 2D
   double const (&momentum_y)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(Equation::MomentumY);
   double const     (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(Equation::Energy);

   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      double const one_radius = 1.0 / (x_block_coordinate + ( (double(i) + 0.5) ) * cell_size);
      unsigned int const index_i = i + CC::FICX();
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         unsigned int const index_j = j + CC::FICY();

         volume_forces[ETI(Equation::Mass)][i][j][0]      -= momentum_x[index_i][index_j][0] * one_radius;
         volume_forces[ETI(Equation::MomentumX)][i][j][0] -= velocity_x[index_i][index_j][0] * momentum_x[index_i][index_j][0] * one_radius;
         volume_forces[ETI(Equation::MomentumY)][i][j][0] -= velocity_x[index_i][index_j][0] * momentum_y[index_i][index_j][0] * one_radius;
         volume_forces[ETI(Equation::Energy)][i][j][0]    -= velocity_x[index_i][index_j][0] * (energy[index_i][index_j][0] + pressure[index_i][index_j][0]) * one_radius;
      }
   }
}
