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
#ifndef WENO5Z_H
#define WENO5Z_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute fluxes according to \cite borges2008
 */
class WENO5Z : public Stencil<WENO5Z> {

   friend Stencil;

   static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

   // Coefficients for WENO5Z scheme
   static constexpr double coef_smoothness_1_  = 13.0/12.0;
   static constexpr double coef_smoothness_2_  = 0.25;

   static constexpr double coef_smoothness_11_ =  1.0;
   static constexpr double coef_smoothness_12_ = -2.0;
   static constexpr double coef_smoothness_13_ =  1.0;
   static constexpr double coef_smoothness_14_ =  1.0;
   static constexpr double coef_smoothness_15_ = -4.0;
   static constexpr double coef_smoothness_16_ =  3.0;

   static constexpr double coef_smoothness_21_ =  1.0;
   static constexpr double coef_smoothness_22_ = -2.0;
   static constexpr double coef_smoothness_23_ =  1.0;
   static constexpr double coef_smoothness_24_ =  1.0;
   static constexpr double coef_smoothness_25_ = -1.0;

   static constexpr double coef_smoothness_31_ =  1.0;
   static constexpr double coef_smoothness_32_ = -2.0;
   static constexpr double coef_smoothness_33_ =  1.0;
   static constexpr double coef_smoothness_34_ =  3.0;
   static constexpr double coef_smoothness_35_ = -4.0;
   static constexpr double coef_smoothness_36_ =  1.0;

   static constexpr double coef_weights_1_ = 0.1;
   static constexpr double coef_weights_2_ = 0.6;
   static constexpr double coef_weights_3_ = 0.3;

   static constexpr double coef_stencils_1_ =  2.0/6.0;
   static constexpr double coef_stencils_2_ = -7.0/6.0;
   static constexpr double coef_stencils_3_ = 11.0/6.0;
   static constexpr double coef_stencils_4_ = -1.0/6.0;
   static constexpr double coef_stencils_5_ =  5.0/6.0;
   static constexpr double coef_stencils_6_ =  2.0/6.0;
   static constexpr double coef_stencils_7_ =  2.0/6.0;
   static constexpr double coef_stencils_8_ =  5.0/6.0;
   static constexpr double coef_stencils_9_ = -1.0/6.0;

   // Number of cells required for upwind and downwind stencils, as well as number of cells downstream of the cell
   static constexpr unsigned int stencil_size_            = 6;
   static constexpr unsigned int downstream_stencil_size_ = 2;

   double ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, const double cell_size) const;

public:
   explicit WENO5Z() = default;
   ~WENO5Z() = default;
   WENO5Z( WENO5Z const& ) = delete;
   WENO5Z& operator=( WENO5Z const& ) = delete;
   WENO5Z( WENO5Z&& ) = delete;
   WENO5Z& operator=( WENO5Z&& ) = delete;
};

#endif // STENCIL_WENO5Z_H
