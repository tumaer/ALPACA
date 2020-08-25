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
#ifndef WENO7_H
#define WENO7_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute fluxes according to \cite Balsara2000.
 */
class WENO7 : public Stencil<WENO7> {

   friend Stencil;

   static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

   // Coefficients for WENO7 scheme
   static constexpr double coef_smoothness_0_01_ = 547.0;
   static constexpr double coef_smoothness_0_02_ = -3882.0;
   static constexpr double coef_smoothness_0_03_ = 4642.0;
   static constexpr double coef_smoothness_0_04_ = -1854.0;
   static constexpr double coef_smoothness_0_06_ = 7043.0;
   static constexpr double coef_smoothness_0_07_ = -17246.0;
   static constexpr double coef_smoothness_0_08_ = 7042.0;
   static constexpr double coef_smoothness_0_10_ = 11003.0;
   static constexpr double coef_smoothness_0_11_ = -9402.0;
   static constexpr double coef_smoothness_0_13_ = 2107.0;

   static constexpr double coef_smoothness_1_01_ = 267.0;
   static constexpr double coef_smoothness_1_02_ = -1642.0;
   static constexpr double coef_smoothness_1_03_ = 1602.0;
   static constexpr double coef_smoothness_1_04_ = -494.0;
   static constexpr double coef_smoothness_1_06_ = 2843.0;
   static constexpr double coef_smoothness_1_07_ = -5966.0;
   static constexpr double coef_smoothness_1_08_ = 1922.0;
   static constexpr double coef_smoothness_1_10_ = 3443.0;
   static constexpr double coef_smoothness_1_11_ = -2522.0;
   static constexpr double coef_smoothness_1_13_ = 547.0;

   static constexpr double coef_smoothness_2_01_ = 547.0;
   static constexpr double coef_smoothness_2_02_ = -2522.0;
   static constexpr double coef_smoothness_2_03_ = 1922.0;
   static constexpr double coef_smoothness_2_04_ = -494.0;
   static constexpr double coef_smoothness_2_06_ = 3443.0;
   static constexpr double coef_smoothness_2_07_ = -5966.0;
   static constexpr double coef_smoothness_2_08_ = 1602.0;
   static constexpr double coef_smoothness_2_10_ = 2843.0;
   static constexpr double coef_smoothness_2_11_ = -1642.0;
   static constexpr double coef_smoothness_2_13_ = 267.0;

   static constexpr double coef_smoothness_3_01_ = 2107.0;
   static constexpr double coef_smoothness_3_02_ = -9402.0;
   static constexpr double coef_smoothness_3_03_ = 7042.0;
   static constexpr double coef_smoothness_3_04_ = -1854.0;
   static constexpr double coef_smoothness_3_06_ = 11003.0;
   static constexpr double coef_smoothness_3_07_ = -17246.0;
   static constexpr double coef_smoothness_3_08_ = 4642.0;
   static constexpr double coef_smoothness_3_10_ = 7043.0;
   static constexpr double coef_smoothness_3_11_ = -3882.0;
   static constexpr double coef_smoothness_3_13_ = 547.0;

   static constexpr double coef_weights_1_ = 1.0 / 35.0;
   static constexpr double coef_weights_2_ = 12.0 / 35.0;
   static constexpr double coef_weights_3_ = 18.0 / 35.0;
   static constexpr double coef_weights_4_ = 4.0 / 35.0;

   static constexpr double coef_stencils_1_ = -3.0 / 12.0;
   static constexpr double coef_stencils_2_ = 13.0 / 12.0;
   static constexpr double coef_stencils_3_ = -23.0 / 12.0;
   static constexpr double coef_stencils_4_ = 25.0 / 12.0;

   static constexpr double coef_stencils_6_ = 1.0 / 12.0;
   static constexpr double coef_stencils_7_ = -5.0 / 12.0;
   static constexpr double coef_stencils_8_ = 13.0 / 12.0;
   static constexpr double coef_stencils_9_ = 3.0 / 12.0;

   static constexpr double coef_stencils_11_ = -1.0 / 12.0;
   static constexpr double coef_stencils_12_ = 7.0 / 12.0;
   static constexpr double coef_stencils_13_ = 7.0 / 12.0;
   static constexpr double coef_stencils_14_ = -1.0 / 12.0;

   static constexpr double coef_stencils_16_ = 3.0 / 12.0;
   static constexpr double coef_stencils_17_ = 13.0 / 12.0;
   static constexpr double coef_stencils_18_ = -5.0 / 12.0;
   static constexpr double coef_stencils_19_ = 1.0 / 12.0;

   // Small values to avoid division by 0, but also to adjust dissipation
   static constexpr double epsilon_weno7_ = 1.0e-10;

   // Number of cells required for upwind and downwind stencils, as well as number of cells downstream of the cell
   static constexpr unsigned int stencil_size_            = 8;
   static constexpr unsigned int downstream_stencil_size_ = 3;

   /**
    * @brief Evaluates the stencil according to a WENO-7 scheme. Also See base class.
    * @note Hotpath function.
    */
   constexpr double ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, double const ) const {
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

      double const s1 = v1 * s11 + v2 * s12 + v3 * s13 + v4 * s14;

      double const s21 = coef_smoothness_1_01_ * v2 + coef_smoothness_1_02_ * v3 + coef_smoothness_1_03_ * v4 + coef_smoothness_1_04_ * v5;
      double const s22 = coef_smoothness_1_06_ * v3 + coef_smoothness_1_07_ * v4 + coef_smoothness_1_08_ * v5;
      double const s23 = coef_smoothness_1_10_ * v4 + coef_smoothness_1_11_ * v5;
      double const s24 = coef_smoothness_1_13_ * v5;

      double const s2 = v2 * s21 + v3 * s22 + v4 * s23 + v5 * s24;

      double const s31 = coef_smoothness_2_01_ * v3 + coef_smoothness_2_02_ * v4 + coef_smoothness_2_03_ * v5 + coef_smoothness_2_04_ * v6;
      double const s32 = coef_smoothness_2_06_ * v4 + coef_smoothness_2_07_ * v5 + coef_smoothness_2_08_ * v6;
      double const s33 = coef_smoothness_2_10_ * v5 + coef_smoothness_2_11_ * v6;
      double const s34 = coef_smoothness_2_13_ * v6;

      double const s3 = v3 * s31 + v4 * s32 + v5 * s33 + v6 * s34;

      double const s41 = coef_smoothness_3_01_ * v4 + coef_smoothness_3_02_ * v5 + coef_smoothness_3_03_ * v6 + coef_smoothness_3_04_ * v7;
      double const s42 = coef_smoothness_3_06_ * v5 + coef_smoothness_3_07_ * v6 + coef_smoothness_3_08_ * v7;
      double const s43 = coef_smoothness_3_10_ * v6 + coef_smoothness_3_11_ * v7;
      double const s44 = coef_smoothness_3_13_ * v7;

      double const s4 = v4 * s41 + v5 * s42 + v6 * s43 + v7 * s44;

      // Compute weights
      // NOTE: The epsilon value is used here explicitly to avoid compiler optimizations when the epsilon is added directly to s_i.
      //       This could lead to undesired behavior in case the values s1, ..., s4 are of similar magnitude.
      //       Then, it cannot guaranteed anymore that a division by zero is avoided.
      double const a1 = coef_weights_1_ / ( ( s1 + epsilon_weno7_ ) * ( s1 + epsilon_weno7_ ) );
      double const a2 = coef_weights_2_ / ( ( s2 + epsilon_weno7_ ) * ( s2 + epsilon_weno7_ ) );
      double const a3 = coef_weights_3_ / ( ( s3 + epsilon_weno7_ ) * ( s3 + epsilon_weno7_ ) );
      double const a4 = coef_weights_4_ / ( ( s4 + epsilon_weno7_ ) * ( s4 + epsilon_weno7_ ) );

      double const one_a_sum = 1.0 / ( a1 + a2 + a3 + a4 );

      double const w1 = a1 * one_a_sum;
      double const w2 = a2 * one_a_sum;
      double const w3 = a3 * one_a_sum;
      double const w4 = a4 * one_a_sum;

      // Return weighted average
      return w1 * ( coef_stencils_1_ * v1 + coef_stencils_2_ * v2 + coef_stencils_3_ * v3 + coef_stencils_4_ * v4 ) + w2 * ( coef_stencils_6_ * v2 + coef_stencils_7_ * v3 + coef_stencils_8_ * v4 + coef_stencils_9_ * v5 ) + w3 * ( coef_stencils_11_ * v3 + coef_stencils_12_ * v4 + coef_stencils_13_ * v5 + coef_stencils_14_ * v6 ) + w4 * ( coef_stencils_16_ * v4 + coef_stencils_17_ * v5 + coef_stencils_18_ * v6 + coef_stencils_19_ * v7 );
   }

public:
   explicit constexpr WENO7() = default;
   ~WENO7()                   = default;
   WENO7( WENO7 const& )      = delete;
   WENO7& operator=( WENO7 const& ) = delete;
   WENO7( WENO7&& )                 = delete;
   WENO7& operator=( WENO7&& ) = delete;
};

#endif// STENCIL_WENO7_H
