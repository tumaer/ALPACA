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
#ifndef WENOCU6_H
#define WENOCU6_H

#include "stencils/stencil.h"

/**
 * @brief Discretization of the SpatialReconstructionStencil class to compute fluxes according to \cite Hu2010.
 */
class WENOCU6 : public Stencil<WENOCU6> {

   friend Stencil;

   static constexpr StencilType stencil_type_ = StencilType::Reconstruction;

   // Coefficients for WENOCU-6 scheme
   static constexpr double coef_smoothness_1_  = 13.0;
   static constexpr double coef_smoothness_2_  = 3.0;

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

   static constexpr double coef_smoothness_411_ =   341.0/5760.0;
   static constexpr double coef_smoothness_412_ = -2785.0/5760.0;
   static constexpr double coef_smoothness_413_ = -2590.0/5760.0;
   static constexpr double coef_smoothness_414_ =  6670.0/5760.0;
   static constexpr double coef_smoothness_415_ = -1895.0/5760.0;
   static constexpr double coef_smoothness_416_ =   259.0/5760.0;

   static constexpr double coef_smoothness_421_ =  -1.0/16.0;
   static constexpr double coef_smoothness_422_ =  12.0/16.0;
   static constexpr double coef_smoothness_423_ = -22.0/16.0;
   static constexpr double coef_smoothness_424_ =  12.0/16.0;
   static constexpr double coef_smoothness_425_ =  -1.0/16.0;

   static constexpr double coef_smoothness_431_ =  -5.0/144.0;
   static constexpr double coef_smoothness_432_ = -11.0/144.0;
   static constexpr double coef_smoothness_433_ =  70.0/144.0;
   static constexpr double coef_smoothness_434_ = -94.0/144.0;
   static constexpr double coef_smoothness_435_ =  47.0/144.0;
   static constexpr double coef_smoothness_436_ =  -7.0/144.0;

   static constexpr double coef_smoothness_441_ =  1.0/24.0;
   static constexpr double coef_smoothness_442_ = -4.0/24.0;
   static constexpr double coef_smoothness_443_ =  6.0/24.0;
   static constexpr double coef_smoothness_444_ = -4.0/24.0;
   static constexpr double coef_smoothness_445_ =  1.0/24.0;

   static constexpr double coef_smoothness_461_ =  -1.0/120.0;
   static constexpr double coef_smoothness_462_ =   5.0/120.0;
   static constexpr double coef_smoothness_463_ = -10.0/120.0;
   static constexpr double coef_smoothness_464_ =  10.0/120.0;
   static constexpr double coef_smoothness_465_ =  -5.0/120.0;
   static constexpr double coef_smoothness_466_ =   1.0/120.0;

   static constexpr double coef_smoothness_41_ =         12.0;
   static constexpr double coef_smoothness_42_ =         52.0;
   static constexpr double coef_smoothness_43_ =          6.0;
   static constexpr double coef_smoothness_44_ =      37548.0/80.0;
   static constexpr double coef_smoothness_45_ =        252.0/5.0;
   static constexpr double coef_smoothness_46_ =         12.0/8.0;
   static constexpr double coef_smoothness_47_ =    1051404.0/140.0;
   static constexpr double coef_smoothness_48_ =     169524.0/224.0;
   static constexpr double coef_smoothness_49_ = 3028045620.0/16128.0;

   static constexpr double coef_smoothness_511_ = 1.0/6.0;
   static constexpr double coef_smoothness_512_ = 4.0/6.0;
   static constexpr double coef_smoothness_513_ = 1.0/6.0;

   static constexpr double weight_const_ = 20.0;

   static constexpr double coef_weights_1_ = 0.05;
   static constexpr double coef_weights_2_ = 0.45;
   static constexpr double coef_weights_3_ = 0.45;
   static constexpr double coef_weights_4_ = 0.05;

   static constexpr double coef_stencils_01_ =  2.0/6.0;
   static constexpr double coef_stencils_02_ = -7.0/6.0;
   static constexpr double coef_stencils_03_ = 11.0/6.0;
   static constexpr double coef_stencils_04_ = -1.0/6.0;
   static constexpr double coef_stencils_05_ =  5.0/6.0;
   static constexpr double coef_stencils_06_ =  2.0/6.0;
   static constexpr double coef_stencils_07_ =  2.0/6.0;
   static constexpr double coef_stencils_08_ =  5.0/6.0;
   static constexpr double coef_stencils_09_ = -1.0/6.0;
   static constexpr double coef_stencils_10_ = 11.0/6.0;
   static constexpr double coef_stencils_11_ = -7.0/6.0;
   static constexpr double coef_stencils_12_ =  2.0/6.0;

   // Number of cells required for upwind and downwind stencils, as well as number of cells downstream of the cell
   static constexpr unsigned int stencil_size_            = 6;
   static constexpr unsigned int downstream_stencil_size_ = 2;

   double ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, double const cell_size) const;

public:
   explicit WENOCU6() = default;
   ~WENOCU6() = default;
   WENOCU6( WENOCU6 const& ) = delete;
   WENOCU6& operator=( WENOCU6 const& ) = delete;
   WENOCU6( WENOCU6&& ) = delete;
   WENOCU6& operator=( WENOCU6&& ) = delete;
};

#endif // STENCIL_WENOCU6_H
