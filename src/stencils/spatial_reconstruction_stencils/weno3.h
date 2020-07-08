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
#ifndef WENO3_H
#define WENO3_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute fluxes according to \cite Shu1999.
 */
class WENO3 : public Stencil<WENO3> {

   friend Stencil;

   static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

   // Coefficients for WENO3 scheme
   static constexpr double coef_smoothness_11_ = -1.0;
   static constexpr double coef_smoothness_12_ =  1.0;

   static constexpr double coef_smoothness_21_ = -1.0;
   static constexpr double coef_smoothness_22_ =  1.0;

   static constexpr double coef_weights_1_ = 1.0/3.0;
   static constexpr double coef_weights_2_ = 2.0/3.0;

   static constexpr double coef_stencils_1_ = -1.0/2.0;
   static constexpr double coef_stencils_2_ =  3.0/2.0;
   static constexpr double coef_stencils_3_ =  1.0/2.0;
   static constexpr double coef_stencils_4_ =  1.0/2.0;

   // Small values to avoid division by 0, but also to adjust dissipation. Optimized according to F. Schranner (same as WENO5)
   static constexpr double epsilon_1_ = 1.0e-6;
   static constexpr double epsilon_2_ = 1.0e-15;

   // Number of cells required for upwind and downwind stencils, as well as number of cells downstream of the cell
   static constexpr unsigned int stencil_size_            = 4;
   static constexpr unsigned int downstream_stencil_size_ = 1;

   /**
    * @brief Evaluates the stencil according to a WENO-3 scheme. Also See base class.
    * @note Hotpath function.
    */
   constexpr double ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, double const ) const {
      // Assign values to v_i to make it easier to read
      double const v1 = array[downstream_stencil_size_ + evaluation_properties[0] - 1 * evaluation_properties[1]];
      double const v2 = array[downstream_stencil_size_ + evaluation_properties[0]];
      double const v3 = array[downstream_stencil_size_ + evaluation_properties[0] + 1 * evaluation_properties[1]];

      // Compute smoothness indicators s_i
      double const s11 = coef_smoothness_11_ * v1 + coef_smoothness_12_ * v2;
      double const s1  = s11 * s11 + epsilon_;

      double const s21 = coef_smoothness_21_ * v2 + coef_smoothness_22_ * v3;
      double const s2  = s21 * s21 + epsilon_;

      // Compute weights
      double const a1 = coef_weights_1_ / (s1 * s1);
      double const a2 = coef_weights_2_ / (s2 * s2);

      double const one_a_sum = 1.0 / (a1 + a2);

      double const w1 = a1 * one_a_sum;
      double const w2 = a2 * one_a_sum;

      // Return weighted average
      return  w1 * (coef_stencils_1_ * v1 + coef_stencils_2_ * v2)
            + w2 * (coef_stencils_3_ * v2 + coef_stencils_4_ * v3);
   }

public:
   explicit constexpr WENO3() = default;
   ~WENO3() = default;
   WENO3( WENO3 const& ) = delete;
   WENO3& operator=( WENO3 const& ) = delete;
   WENO3( WENO3&& ) = delete;
   WENO3& operator=( WENO3&& ) = delete;
};

#endif // STENCIL_WENO3_H
