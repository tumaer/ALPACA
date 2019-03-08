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
#include "weno-ao53.h"

#include <stdexcept>
#include <cmath>

/**
 * @brief Computes the flux at one cell face according to WENO-AO(5,3) adaptive-order scheme from Balsara (2016). Also See base class.
 */
double WENOAO53::ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, const double cell_size) const {

#ifndef PERFORMANCE
   (void)cell_size;

   // Output error in case something went wrong with the stencil size
   if(array.size() < stencil_size_) {
      throw std::logic_error("Stencil size in WENO-AO53 is longer than provided array");
   }
#endif

   // Assign values to v_i to make it easier to read
   const double v1 = array[downstream_stencil_size_ + evaluation_properties[0] - 2 * evaluation_properties[1]];
   const double v2 = array[downstream_stencil_size_ + evaluation_properties[0] - 1 * evaluation_properties[1]];
   const double v3 = array[downstream_stencil_size_ + evaluation_properties[0]];
   const double v4 = array[downstream_stencil_size_ + evaluation_properties[0] + 1 * evaluation_properties[1]];
   const double v5 = array[downstream_stencil_size_ + evaluation_properties[0] + 2 * evaluation_properties[1]];

   // Compute Legendre coefficents u according to Eq. (2.3-2.5), Eq. (2.16) amd smoothness indicators beta according to Eq. (2.6) and Eq. (2.19)
   const double u_r3_11 = coef_smoothness_1_1_ * v1 + coef_smoothness_1_2_ * v2 + coef_smoothness_1_3_ * v3;
   const double u_r3_12 = coef_smoothness_1_4_ * v1 + coef_smoothness_1_5_ * v2 + coef_smoothness_1_6_ * v3;
   const double beta_r3_1 = coef_smoothness_1_ * u_r3_11 * u_r3_11 + coef_smoothness_2_ * u_r3_12 * u_r3_12;

   const double u_r3_21 = coef_smoothness_2_1_ * v2 + coef_smoothness_2_2_ * v4;
   const double u_r3_22 = coef_smoothness_2_3_ * v2 + coef_smoothness_2_4_ * v3 + coef_smoothness_2_5_ * v4;
   const double beta_r3_2 = coef_smoothness_1_ * u_r3_21 * u_r3_21 + coef_smoothness_2_ * u_r3_22 * u_r3_22;

   const double u_r3_31 = coef_smoothness_3_1_ * v3 + coef_smoothness_3_2_ * v4 + coef_smoothness_3_3_ * v5;
   const double u_r3_32 = coef_smoothness_3_4_ * v3 + coef_smoothness_3_5_ * v4 + coef_smoothness_3_6_ * v5;
   const double beta_r3_3 = coef_smoothness_1_ * u_r3_31 * u_r3_31 + coef_smoothness_2_ * u_r3_32 * u_r3_32;

   const double u_r5_31 = coef_smoothness_5_01_ * v1 + coef_smoothness_5_02_ * v2 + coef_smoothness_5_03_ * v4 + coef_smoothness_5_04_ * v5;
   const double u_r5_32 = coef_smoothness_5_05_ * v1 + coef_smoothness_5_06_ * v2 + coef_smoothness_5_07_ * v3 + coef_smoothness_5_08_ * v4 + coef_smoothness_5_09_ * v5;
   const double u_r5_33 = coef_smoothness_5_10_ * v1 + coef_smoothness_5_11_ * v2 + coef_smoothness_5_12_ * v4 + coef_smoothness_5_13_ * v5;
   const double u_r5_34 = coef_smoothness_5_14_ * v1 + coef_smoothness_5_15_ * v2 + coef_smoothness_5_16_ * v3 + coef_smoothness_5_17_ * v4 + coef_smoothness_5_18_ * v5;

   const double beta_r5_3 = coef_smoothness_1_ * (u_r5_31 + coef_smoothness_weight_5_1_ * u_r5_33) * (u_r5_31 + coef_smoothness_weight_5_1_ * u_r5_33)
      + coef_smoothness_2_ * (u_r5_32 + coef_smoothness_weight_5_2_ * u_r5_34) * (u_r5_32 + coef_smoothness_weight_5_2_ * u_r5_34)
      + coef_smoothness_3_ * u_r5_33 * u_r5_33 + coef_smoothness_4_ * u_r5_34 * u_r5_34;

   // Compute normalized weights
   // Eq. (3.6)
   const double tau = (std::abs( beta_r5_3 - beta_r3_1 ) + std::abs( beta_r5_3 - beta_r3_2 ) + std::abs( beta_r5_3 - beta_r3_3 )) * one_third_;

   // Eq. (3.7a) Note: Balsara et al. suggest an epsilon value of 1e-12 to minimize the influence. We use machine precision instead.
   const double a1 = linear_weight_r3_1_ * ( 1.0 + (tau * tau) / ((beta_r3_1 + epsilon_) * (beta_r3_1 + epsilon_)) );
   const double a2 = linear_weight_r3_2_ * ( 1.0 + (tau * tau) / ((beta_r3_2 + epsilon_) * (beta_r3_2 + epsilon_)) );
   const double a3 = linear_weight_r3_3_ * ( 1.0 + (tau * tau) / ((beta_r3_3 + epsilon_) * (beta_r3_3 + epsilon_)) );
   const double a5 = linear_weight_r5_3_ * ( 1.0 + (tau * tau) / ((beta_r5_3 + epsilon_) * (beta_r5_3 + epsilon_)) );

   const double one_a_sum = 1.0 / (a1 + a2 + a3 + a5);

   const double w1 = a1 * one_a_sum;
   const double w2 = a2 * one_a_sum;
   const double w3 = a3 * one_a_sum;
   const double w5 = a5 * one_a_sum;

   // Compute coefficients of the Legendre basis polynomial according to Eq. (3.10)
   const double tmp1 = w5 * one_over_linear_weight_r5_3_;

   const double u0 = v3;
   const double u1 = tmp1 * (u_r5_31 - linear_weight_r3_1_ * u_r3_11 - linear_weight_r3_2_ * u_r3_21 - linear_weight_r3_3_ * u_r3_31) + w1 * u_r3_11 + w2 * u_r3_21 + w3 * u_r3_31;
   const double u2 = tmp1 * (u_r5_32 - linear_weight_r3_1_ * u_r3_12 - linear_weight_r3_2_ * u_r3_22 - linear_weight_r3_3_ * u_r3_32) + w1 * u_r3_12 + w2 * u_r3_22 + w3 * u_r3_32;
   const double u3 = tmp1 * u_r5_33;
   const double u4 = tmp1 * u_r5_34;

   // Return value of reconstructed polynomial according to Eq. (3.11) Note: The polynomial consturcted in the paper is always evaluated for x = 0.5 (cell face).
   return  u0 + legendre_1_ * u1 + legendre_2_ * u2 + legendre_3_ * u3 + legendre_4_ * u4;
}
