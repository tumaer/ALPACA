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
#include "weno7.h"

#include <stdexcept>

/**
 * @brief Computes the flux at one cell face according to used WENO-7 scheme. Also See base class.
 * @note Hotpath function.
 */
double WENO7::ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, double const cell_size) const {

#ifndef PERFORMANCE
   (void)cell_size;

   // Output error in case something went wrong with the stencil size
   if(array.size() < stencil_size_) {
      throw std::logic_error("Stencil size in WENO7 is longer than provided Array");
   }
#endif

   // Assign values to v_i to make it easier to read
   double const v1 = array[downstream_stencil_size_ + evaluation_properties[0] - 3 * evaluation_properties[1]];
   double const v2 = array[downstream_stencil_size_ + evaluation_properties[0] - 2 * evaluation_properties[1]];
   double const v3 = array[downstream_stencil_size_ + evaluation_properties[0] - 1 * evaluation_properties[1]];
   double const v4 = array[downstream_stencil_size_ + evaluation_properties[0]];
   double const v5 = array[downstream_stencil_size_ + evaluation_properties[0] + 1 * evaluation_properties[1]];
   double const v6 = array[downstream_stencil_size_ + evaluation_properties[0] + 2 * evaluation_properties[1]];
   double const v7 = array[downstream_stencil_size_ + evaluation_properties[0] + 3 * evaluation_properties[1]];

   // Compute smoothness indicators s_i
   double const s11 = coef_smoothness_0_01_ * v1 + coef_smoothness_0_02_ * v2 + coef_smoothness_0_03_ * v3 + coef_smoothness_0_04_ * v4;
   double const s12 = coef_smoothness_0_06_ * v2 + coef_smoothness_0_07_ * v3 + coef_smoothness_0_08_ * v4;
   double const s13 = coef_smoothness_0_10_ * v3 + coef_smoothness_0_11_ * v4;
   double const s14 = coef_smoothness_0_13_ * v4;

   double const s1 = (v1*s11 + v2*s12 + v3*s13 + v4*s14) + epsilon_weno7_;

   double const s21 = coef_smoothness_1_01_ * v2 + coef_smoothness_1_02_ * v3 + coef_smoothness_1_03_ * v4 + coef_smoothness_1_04_ * v5;
   double const s22 = coef_smoothness_1_06_ * v3 + coef_smoothness_1_07_ * v4 + coef_smoothness_1_08_ * v5;
   double const s23 = coef_smoothness_1_10_ * v4 + coef_smoothness_1_11_ * v5;
   double const s24 = coef_smoothness_1_13_ * v5;

   double const s2 = (v2*s21 + v3*s22 + v4*s23 + v5*s24) + epsilon_weno7_;

   double const s31 = coef_smoothness_2_01_ * v3 + coef_smoothness_2_02_ * v4 + coef_smoothness_2_03_ * v5 + coef_smoothness_2_04_ * v6;
   double const s32 = coef_smoothness_2_06_ * v4 + coef_smoothness_2_07_ * v5 + coef_smoothness_2_08_ * v6;
   double const s33 = coef_smoothness_2_10_ * v5 + coef_smoothness_2_11_ * v6;
   double const s34 = coef_smoothness_2_13_ * v6;

   double const s3 = (v3*s31 + v4*s32 + v5*s33 + v6*s34) + epsilon_weno7_;

   double const s41 = coef_smoothness_3_01_ * v4 + coef_smoothness_3_02_ * v5 + coef_smoothness_3_03_ * v6 + coef_smoothness_3_04_ * v7;
   double const s42 = coef_smoothness_3_06_ * v5 + coef_smoothness_3_07_ * v6 + coef_smoothness_3_08_ * v7;
   double const s43 = coef_smoothness_3_10_ * v6 + coef_smoothness_3_11_ * v7;
   double const s44 = coef_smoothness_3_13_ * v7;

   double const s4 = (v4*s41 + v5*s42 + v6*s43 + v7*s44) + epsilon_weno7_;

   // Compute weights
   double const a1 = coef_weights_1_ / (s1 * s1);
   double const a2 = coef_weights_2_ / (s2 * s2);
   double const a3 = coef_weights_3_ / (s3 * s3);
   double const a4 = coef_weights_4_ / (s4 * s4);

   double const one_a_sum = 1.0 / (a1 + a2 + a3 + a4);

   double const w1 = a1 * one_a_sum;
   double const w2 = a2 * one_a_sum;
   double const w3 = a3 * one_a_sum;
   double const w4 = a4 * one_a_sum;

   // Return weighted average
   return  w1 * (coef_stencils_1_  * v1 + coef_stencils_2_  * v2 + coef_stencils_3_  * v3 + coef_stencils_4_  * v4 )
         + w2 * (coef_stencils_6_  * v2 + coef_stencils_7_  * v3 + coef_stencils_8_  * v4 + coef_stencils_9_  * v5 )
         + w3 * (coef_stencils_11_ * v3 + coef_stencils_12_ * v4 + coef_stencils_13_ * v5 + coef_stencils_14_ * v6 )
         + w4 * (coef_stencils_16_ * v4 + coef_stencils_17_ * v5 + coef_stencils_18_ * v6 + coef_stencils_19_ * v7 );
}
