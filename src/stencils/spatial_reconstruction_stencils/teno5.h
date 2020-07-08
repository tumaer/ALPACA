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
#ifndef TENO5_H
#define TENO5_H

#include "stencils/stencil.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute fluxes according to \cite Fu2016a.
 */
class TENO5 : public Stencil<TENO5> {

   friend Stencil;

   static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

   // Ideal weights and threshold for stencil smoothness of TENO5
   static constexpr double d1_ = 0.6;
   static constexpr double d2_ = 0.3;
   static constexpr double d3_ = 0.1;
   static constexpr double CT_ = 1.0e-5;

   // Coefficients for TENO5 scheme
   static constexpr double coef_smoothness_1_  = 13.0/12.0;
   static constexpr double coef_smoothness_2_  = 0.25;

   static constexpr double coef_smoothness_11_ =  1.0;
   static constexpr double coef_smoothness_12_ = -2.0;
   static constexpr double coef_smoothness_13_ =  1.0;
   static constexpr double coef_smoothness_14_ =  1.0;
   static constexpr double coef_smoothness_15_ = -1.0;

   static constexpr double coef_smoothness_21_ =  1.0;
   static constexpr double coef_smoothness_22_ = -2.0;
   static constexpr double coef_smoothness_23_ =  1.0;
   static constexpr double coef_smoothness_24_ =  3.0;
   static constexpr double coef_smoothness_25_ = -4.0;
   static constexpr double coef_smoothness_26_ =  1.0;

   static constexpr double coef_smoothness_31_ =  1.0;
   static constexpr double coef_smoothness_32_ = -2.0;
   static constexpr double coef_smoothness_33_ =  1.0;
   static constexpr double coef_smoothness_34_ =  1.0;
   static constexpr double coef_smoothness_35_ = -4.0;
   static constexpr double coef_smoothness_36_ =  3.0;

   static constexpr double coef_stencils_1_ = -1.0;
   static constexpr double coef_stencils_2_ =  5.0;
   static constexpr double coef_stencils_3_ =  2.0;
   static constexpr double coef_stencils_4_ =  2.0;
   static constexpr double coef_stencils_5_ =  5.0;
   static constexpr double coef_stencils_6_ = -1.0;
   static constexpr double coef_stencils_7_ =  2.0;
   static constexpr double coef_stencils_8_ = -7.0;
   static constexpr double coef_stencils_9_ = 11.0;

   static constexpr double multiplyer_stencils_ = 1.0/6.0;

   // Number of cells required for upwind and downwind stencils, as well as number of cells downstream of the cell
   static constexpr unsigned int stencil_size_            = 6;
   static constexpr unsigned int downstream_stencil_size_ = 2;

   /**
    * @brief Evaluates the stencil according to a TENO scheme of fifth order. Also See base class.
    * @note Hotpath function.
    */
   constexpr double ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, double const ) const {
      // Assign values to v_i to make it easier to read
      double const v1 = array[downstream_stencil_size_ + evaluation_properties[0] - 2 * evaluation_properties[1]];
      double const v2 = array[downstream_stencil_size_ + evaluation_properties[0] - 1 * evaluation_properties[1]];
      double const v3 = array[downstream_stencil_size_ + evaluation_properties[0]];
      double const v4 = array[downstream_stencil_size_ + evaluation_properties[0] + 1 * evaluation_properties[1]];
      double const v5 = array[downstream_stencil_size_ + evaluation_properties[0] + 2 * evaluation_properties[1]];

      // Compute smoothness indicators si
      double const s11 = coef_smoothness_11_ * v2 + coef_smoothness_12_ * v3 + coef_smoothness_13_ * v4;
      double const s12 = coef_smoothness_14_ * v2 + coef_smoothness_15_ * v4;

      double const s1 = coef_smoothness_1_*s11*s11 + coef_smoothness_2_*s12*s12;

      double const s21 = coef_smoothness_21_ * v3 + coef_smoothness_22_ * v4 + coef_smoothness_23_ * v5;
      double const s22 = coef_smoothness_24_ * v3 + coef_smoothness_25_ * v4 + coef_smoothness_26_ * v5;

      double const s2 = coef_smoothness_1_*s21*s21 + coef_smoothness_2_*s22*s22;

      double const s31 = coef_smoothness_31_ * v1 + coef_smoothness_32_ * v2 + coef_smoothness_33_ * v3;
      double const s32 = coef_smoothness_34_ * v1 + coef_smoothness_35_ * v2 + coef_smoothness_36_ * v3;

      double const s3 = coef_smoothness_1_*s31*s31 + coef_smoothness_2_*s32*s32;

      double const tau5 = Abs( s3 - s2 );

      double a1 = 1.0 + tau5 / (s1 + epsilon_);
      double a2 = 1.0 + tau5 / (s2 + epsilon_);
      double a3 = 1.0 + tau5 / (s3 + epsilon_);

      // NF Calculate a^6 without using std::pow to improve performance drastically
      a1 *= (a1*a1);
      a2 *= (a2*a2);
      a3 *= (a3*a3);
      a1 *= a1;
      a2 *= a2;
      a3 *= a3;

      double const one_a_sum = 1.0 / (a1 + a2 + a3);

      double const b1 = a1 * one_a_sum < CT_ ? 0.0 : 1.0;
      double const b2 = a2 * one_a_sum < CT_ ? 0.0 : 1.0;
      double const b3 = a3 * one_a_sum < CT_ ? 0.0 : 1.0;

      double const Variation1 = coef_stencils_1_ * v2 + coef_stencils_2_ * v3 + coef_stencils_3_ * v4;
      double const Variation2 = coef_stencils_4_ * v3 + coef_stencils_5_ * v4 + coef_stencils_6_ * v5;
      double const Variation3 = coef_stencils_7_ * v1 + coef_stencils_8_ * v2 + coef_stencils_9_ * v3;

      double const w1 = d1_ * b1;
      double const w2 = d2_ * b2;
      double const w3 = d3_ * b3;

      double const one_w_sum = 1.0 / (w1 + w2 + w3);

      double const w1_normalized = w1 * one_w_sum;
      double const w2_normalized = w2 * one_w_sum;
      double const w3_normalized = w3 * one_w_sum;

      return (w1_normalized * Variation1 + w2_normalized * Variation2 + w3_normalized * Variation3) * multiplyer_stencils_;
   }

public:
   explicit constexpr TENO5() = default;
   ~TENO5() = default;
   TENO5( TENO5 const& ) = delete;
   TENO5& operator=( TENO5 const& ) = delete;
   TENO5( TENO5&& ) = delete;
   TENO5& operator=( TENO5&& ) = delete;
};

#endif // STENCIL_TENO5_H
