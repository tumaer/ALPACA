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
#include "weno3.h"

#include <stdexcept>

/**
 * @brief Computes the flux at one cell face according to used WENO-3 scheme. Also See base class.
 */
double WENO3::ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, const double cell_size) const {

#ifndef PERFORMANCE
   (void)cell_size;

   // Output error in case something went wrong with the stencil size
   if(array.size() < stencil_size_) {
      throw std::logic_error("Stencil size in WENO3 is longer than provided Array");
   }
#endif

   // Assign values to v_i to make it easier to read
   const double v1 = array[downstream_stencil_size_ + evaluation_properties[0] - 1 * evaluation_properties[1]];
   const double v2 = array[downstream_stencil_size_ + evaluation_properties[0]];
   const double v3 = array[downstream_stencil_size_ + evaluation_properties[0] + 1 * evaluation_properties[1]];

   // Compute smoothness indicators s_i
   const double s11 = coef_smoothness_11_ * v1 + coef_smoothness_12_ * v2;
   const double s1  = s11 * s11 + epsilon_;

   const double s21 = coef_smoothness_21_ * v2 + coef_smoothness_22_ * v3;
   const double s2  = s21 * s21 + epsilon_;

   // Compute weights
   const double a1 = coef_weights_1_ / (s1 * s1);
   const double a2 = coef_weights_2_ / (s2 * s2);

   const double one_a_sum = 1.0 / (a1 + a2);

   const double w1 = a1 * one_a_sum;
   const double w2 = a2 * one_a_sum;

   // Return weighted average
   return  w1 * (coef_stencils_1_ * v1 + coef_stencils_2_ * v2)
      + w2 * (coef_stencils_3_ * v2 + coef_stencils_4_ * v3);
}
