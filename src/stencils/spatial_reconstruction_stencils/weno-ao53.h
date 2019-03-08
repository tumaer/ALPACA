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
#ifndef WENOAO53_H
#define WENOAO53_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute fluxes according to \cite Balsara2016.
 */
class WENOAO53 : public Stencil<WENOAO53> {

   friend Stencil;

   static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

   // Parameters to control the order reduction (recommended range [0.85,0.95])
   static constexpr double gamma_hi_ = 0.85;
   static constexpr double gamma_lo_ = 0.85;

   // Coefficients for WENOAO53 scheme
   static constexpr double coef_smoothness_1_  = 1.0;
   static constexpr double coef_smoothness_2_  = 13.0/3.0;
   static constexpr double coef_smoothness_3_  = 781.0/20.0;
   static constexpr double coef_smoothness_4_  = 1421461.0/2275.0;

   static constexpr double coef_smoothness_1_1_ =  1.0/2.0;
   static constexpr double coef_smoothness_1_2_ = -2.0;
   static constexpr double coef_smoothness_1_3_ =  3.0/2.0;
   static constexpr double coef_smoothness_1_4_ =  1.0/2.0;
   static constexpr double coef_smoothness_1_5_ = -1.0;
   static constexpr double coef_smoothness_1_6_ =  1.0/2.0;

   static constexpr double coef_smoothness_2_1_ = -1.0/2.0;
   static constexpr double coef_smoothness_2_2_ =  1.0/2.0;
   static constexpr double coef_smoothness_2_3_ =  1.0/2.0;
   static constexpr double coef_smoothness_2_4_ = -1.0;
   static constexpr double coef_smoothness_2_5_ =  1.0/2.0;

   static constexpr double coef_smoothness_3_1_ = -3.0/2.0;
   static constexpr double coef_smoothness_3_2_ =  2.0;
   static constexpr double coef_smoothness_3_3_ = -1.0/2.0;
   static constexpr double coef_smoothness_3_4_ =  1.0/2.0;
   static constexpr double coef_smoothness_3_5_ = -1.0;
   static constexpr double coef_smoothness_3_6_ =  1.0/2.0;

   static constexpr double coef_smoothness_5_01_ =  11.0/120.0;
   static constexpr double coef_smoothness_5_02_ = -82.0/120.0;
   static constexpr double coef_smoothness_5_03_ =  82.0/120.0;
   static constexpr double coef_smoothness_5_04_ = -11.0/120.0;
   static constexpr double coef_smoothness_5_05_ = -3.0/56.0;
   static constexpr double coef_smoothness_5_06_ =  40.0/56.0;
   static constexpr double coef_smoothness_5_07_ = -74.0/56.0;
   static constexpr double coef_smoothness_5_08_ =  40.0/56.0;
   static constexpr double coef_smoothness_5_09_ = -3.0/56.0;
   static constexpr double coef_smoothness_5_10_ = -1.0/12.0;
   static constexpr double coef_smoothness_5_11_ =  2.0/12.0;
   static constexpr double coef_smoothness_5_12_ = -2.0/12.0;
   static constexpr double coef_smoothness_5_13_ =  1.0/12.0;
   static constexpr double coef_smoothness_5_14_ =  1.0/24.0;
   static constexpr double coef_smoothness_5_15_ = -4.0/24.0;
   static constexpr double coef_smoothness_5_16_ =  6.0/24.0;
   static constexpr double coef_smoothness_5_17_ = -4.0/24.0;
   static constexpr double coef_smoothness_5_18_ =  1.0/24.0;

   static constexpr double coef_smoothness_weight_5_1_ =  1.0/10.0;
   static constexpr double coef_smoothness_weight_5_2_ =  123.0/455.0;

   // Linear weights according to Eq. (3.5)
   static constexpr double linear_weight_r3_1_ = (1.0 - gamma_hi_) * (1.0 - gamma_lo_) * 0.5;
   static constexpr double linear_weight_r3_2_ = (1.0 - gamma_hi_) * gamma_lo_;
   static constexpr double linear_weight_r3_3_ = linear_weight_r3_1_;
   static constexpr double linear_weight_r5_3_ = gamma_hi_;
   static constexpr double one_over_linear_weight_r5_3_ = 1.0 / linear_weight_r5_3_;

   // Legendre polynomials evaluated at cell face x = 1/2
   static constexpr double legendre_1_ = 1.0/2.0;
   static constexpr double legendre_2_ = 1.0/6.0;
   static constexpr double legendre_3_ = 1.0/20.0;
   static constexpr double legendre_4_ = 1.0/70.0;

   // Precompiled constant
   static constexpr double one_third_  = 1.0 / 3.0;

   // Number of cells required for upwind and downwind stencils, as well as number of cells downstream of the cell
   static constexpr unsigned int stencil_size_            = 6;
   static constexpr unsigned int downstream_stencil_size_ = 2;

   double ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, const double cell_size) const;

public:
   explicit WENOAO53() = default;
   ~WENOAO53() = default;
   WENOAO53( WENOAO53 const& ) = delete;
   WENOAO53& operator=( WENOAO53 const& ) = delete;
   WENOAO53( WENOAO53&& ) = delete;
   WENOAO53& operator=( WENOAO53&& ) = delete;
};

#endif // STENCIL_WENOAO53_H
