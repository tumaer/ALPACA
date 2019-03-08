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
#ifndef SPACE_FILLING_CURVES_H
#define SPACE_FILLING_CURVES_H

#include <array>

/* Nomenclature according to \cite Bader 2013
 * _x implies "x bar"
 */
enum class HilbertPosition : unsigned short {xyz, x_y_z, yzx, y_z_x, zxy, z_x_y, _xy_z, _x_yz, _yz_x, _y_zx, _zx_y,  _z_xy};

/**
 * @brief The SpaceFillingCurves namespace provides space filling curve Information %Currently only for Hilbert Curve% for the Load Balancing.
 * $PLEASE NOTE: SPACEFILLING CURVE LOAD BALANCING IS ONLY USEFUL IN 3D SIMULATIONS. HENCE THIS CLASS ONLY PROVIDES 3D VARIANTS$
 */
namespace SpaceFillingCurves {

/**
 * @brief Gives the Hilbert Order for the given position.
 * @param position The current position of the cube in the Hilbert Curve.
 * @return The traversal order.
 */
constexpr std::array<unsigned short,8> GetHilbertOrder( const HilbertPosition position ) {
   switch( position ) {
      case HilbertPosition::xyz:   return {0,2,3,1, 5,7,6,4};
      case HilbertPosition::x_y_z: return {6,4,5,7, 3,1,0,2};
      case HilbertPosition::yzx:   return {0,1,5,4, 6,7,3,2};
      case HilbertPosition::y_z_x: return {6,7,3,2, 0,1,5,4};

      case HilbertPosition::zxy:   return {0,4,6,2, 3,7,5,1};
      case HilbertPosition::z_x_y: return {6,2,0,4, 5,1,3,7};
      case HilbertPosition::_xy_z: return {5,7,6,4, 0,2,3,1};
      case HilbertPosition::_x_yz: return {3,1,0,2, 6,4,5,7};

      case HilbertPosition::_yz_x: return {5,4,0,1, 3,2,6,7};
      case HilbertPosition::_y_zx: return {3,2,6,7, 5,4,0,1};
      case HilbertPosition::_zx_y: return {5,1,3,7, 6,2,0,4};
      case HilbertPosition::_z_xy: return {3,7,5,1, 0,4,6,2};
   }
   return {0,1,2,3, 4,5,6,7};
}

/**
 * @brief Gives the Hilbert Replacement for the given Position.
 * @param position The current position of the cube in the Hilbert Curve.
 * @return The replacement order.
 */
constexpr std::array<HilbertPosition,8> GetHilbertReplacement( const HilbertPosition position ) {
   switch( position ) {
      case HilbertPosition::xyz:   return { HilbertPosition::yzx,    HilbertPosition::zxy,    HilbertPosition::zxy,    HilbertPosition::_x_yz,
                                            HilbertPosition::_x_yz,  HilbertPosition::_zx_y,  HilbertPosition::_zx_y,  HilbertPosition::y_z_x};

      case HilbertPosition::x_y_z: return { HilbertPosition::y_z_x,  HilbertPosition::z_x_y,  HilbertPosition::z_x_y,  HilbertPosition::_xy_z,
                                            HilbertPosition::_xy_z,  HilbertPosition::_z_xy,  HilbertPosition::_z_xy,  HilbertPosition::yzx  };

      case HilbertPosition::yzx:   return { HilbertPosition::zxy,    HilbertPosition::xyz,    HilbertPosition::xyz,    HilbertPosition::_yz_x,
                                            HilbertPosition::_yz_x,  HilbertPosition::x_y_z,  HilbertPosition::x_y_z,  HilbertPosition::_z_xy};

      case HilbertPosition::y_z_x: return { HilbertPosition::z_x_y,  HilbertPosition::x_y_z,  HilbertPosition::x_y_z,  HilbertPosition::_y_zx,
                                            HilbertPosition::_y_zx,  HilbertPosition::xyz,    HilbertPosition::xyz,    HilbertPosition::_zx_y};


      case HilbertPosition::zxy:   return { HilbertPosition::xyz,    HilbertPosition::yzx,    HilbertPosition::yzx,    HilbertPosition::z_x_y,
                                            HilbertPosition::z_x_y,  HilbertPosition::_y_zx,  HilbertPosition::_y_zx,  HilbertPosition::_xy_z};

      case HilbertPosition::z_x_y: return { HilbertPosition::x_y_z,  HilbertPosition::y_z_x,  HilbertPosition::y_z_x,  HilbertPosition::zxy,
                                            HilbertPosition::zxy,    HilbertPosition::_yz_x,  HilbertPosition::_yz_x,  HilbertPosition::_x_yz};

      case HilbertPosition::_xy_z: return { HilbertPosition::_yz_x,  HilbertPosition::_zx_y,  HilbertPosition::_zx_y,  HilbertPosition::x_y_z,
                                            HilbertPosition::x_y_z,  HilbertPosition::zxy,    HilbertPosition::zxy,    HilbertPosition::_y_zx};

      case HilbertPosition::_x_yz: return { HilbertPosition::_y_zx,  HilbertPosition::_z_xy,  HilbertPosition::_z_xy,  HilbertPosition::xyz,
                                            HilbertPosition::xyz,    HilbertPosition::z_x_y,  HilbertPosition::z_x_y,  HilbertPosition::_yz_x};


      case HilbertPosition::_yz_x: return { HilbertPosition::_zx_y,  HilbertPosition::_xy_z,  HilbertPosition::_xy_z,  HilbertPosition::yzx,
                                            HilbertPosition::yzx,    HilbertPosition::_x_yz,  HilbertPosition::_x_yz,  HilbertPosition::z_x_y};

      case HilbertPosition::_y_zx: return { HilbertPosition::_z_xy,  HilbertPosition::_x_yz,  HilbertPosition::_x_yz,  HilbertPosition::y_z_x,
                                            HilbertPosition::y_z_x,  HilbertPosition::_xy_z,  HilbertPosition::_xy_z,  HilbertPosition::zxy  };

      case HilbertPosition::_zx_y: return { HilbertPosition::_xy_z,  HilbertPosition::_yz_x,  HilbertPosition::_yz_x,  HilbertPosition::_z_xy,
                                            HilbertPosition::_z_xy,  HilbertPosition::y_z_x,  HilbertPosition::y_z_x,  HilbertPosition::xyz  };

      case HilbertPosition::_z_xy: return { HilbertPosition::_x_yz,  HilbertPosition::_y_zx,  HilbertPosition::_y_zx,  HilbertPosition::_zx_y,
                                            HilbertPosition::_zx_y,  HilbertPosition::yzx,    HilbertPosition::yzx,    HilbertPosition::x_y_z};
   }

   return {{ HilbertPosition::yzx,    HilbertPosition::zxy,    HilbertPosition::zxy,    HilbertPosition::_x_yz,
              HilbertPosition::_x_yz,  HilbertPosition::_zx_y,  HilbertPosition::_zx_y,  HilbertPosition::y_z_x}};
}
}

#endif // SPACE_FILLING_CURVES_H
