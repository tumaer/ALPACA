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
#ifndef GEOMETRY_CALCULATOR_MARCHING_CUBES_H
#define GEOMETRY_CALCULATOR_MARCHING_CUBES_H

#include "geometry_calculator.h"

/**
 * @brief The GeometryCalculatorMarchingCubes class calculates volume fractions and cell face apertures based on the marching-cube
 *        algorithm as presented in \cite Lauer2012.
 * @note This class is derived and inherits from the abstract class GeometryCalculator.
 */
class GeometryCalculatorMarchingCubes : public GeometryCalculator<GeometryCalculatorMarchingCubes> {

private:
   /**
    * Bool that indicates whether geometric quantities are calculated cell based. The defaults setting is false, and geometric calculations
    * are performed sub-cell based. I.e., A 3-D cell is split into 8 sub cells.
    */
   static constexpr bool cell_based_geometry_calculations_ = false;

public:
   explicit GeometryCalculatorMarchingCubes()                                = default;
   ~GeometryCalculatorMarchingCubes()                                        = default;
   GeometryCalculatorMarchingCubes( GeometryCalculatorMarchingCubes const& ) = delete;
   GeometryCalculatorMarchingCubes& operator=( GeometryCalculatorMarchingCubes const& ) = delete;
   GeometryCalculatorMarchingCubes( GeometryCalculatorMarchingCubes&& )                 = delete;
   GeometryCalculatorMarchingCubes& operator=( GeometryCalculatorMarchingCubes&& ) = delete;

   std::array<double, 6> ComputeCellFaceApertureImplementation( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                                unsigned int const i, unsigned int const j, unsigned int const k, std::int8_t const material_sign = 1 ) const;

   double ComputeVolumeFractionImplementation( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                               unsigned int const i, unsigned int const j, unsigned int const k,
                                               std::int8_t const material_sign = 1 ) const;
};

#endif//GEOMETRY_CALCULATOR_MARCHING_CUBES_H
