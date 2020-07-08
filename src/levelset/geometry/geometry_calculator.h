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
#ifndef GEOMETRY_CALCULATOR_H
#define GEOMETRY_CALCULATOR_H

#include "user_specifications/compile_time_constants.h"
#include <vector>
#include <enums/geometry_settings.h>
#include "topology/node.h"

/**
 * @brief The GeometryCalculator class provides functions for basic geometric calculations (normals, curvature, etc.) based on the
 *        level-set field and serves as abstract class for more sophisticated geometrical functions (apertures, volume fraction, etc.).
 */
template<typename DerivedGeometryCalculator>
class GeometryCalculator {

public:
   GeometryCalculator()                            = default;
   ~GeometryCalculator()                           = default;
   GeometryCalculator( GeometryCalculator const& ) = delete;
   GeometryCalculator& operator=( GeometryCalculator const& ) = delete;
   GeometryCalculator( GeometryCalculator&& )                 = delete;
   GeometryCalculator& operator=( GeometryCalculator&& ) = delete;

   std::array<double, 6> ComputeCellFaceAperture( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k, std::int8_t const material_sign = 1 ) const {
      return static_cast<const DerivedGeometryCalculator*>( this )->ComputeCellFaceApertureImplementation( levelset, i, j, k, material_sign );
   }

   double ComputeVolumeFraction( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k, std::int8_t const material_sign = 1 ) const {
      return static_cast<const DerivedGeometryCalculator*>( this )->ComputeVolumeFractionImplementation( levelset, i, j, k, material_sign );
   }
};

//box size for subcell and cut-cell computation
static constexpr unsigned int subcell_box_size_x = 3;
static constexpr unsigned int subcell_box_size_y = CC::DIM() != Dimension::One ? 3 : 1;
static constexpr unsigned int subcell_box_size_z = CC::DIM() == Dimension::Three ? 3 : 1;

static constexpr unsigned int cell_box_size_x = 2;
static constexpr unsigned int cell_box_size_y = CC::DIM() != Dimension::One ? 2 : 1;
static constexpr unsigned int cell_box_size_z = CC::DIM() == Dimension::Three ? 2 : 1;

void GetLevelsetAtCellCorners( double ( &cell_corner_levelset )[cell_box_size_x][cell_box_size_y][cell_box_size_z],
                               double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()],
                               unsigned int const i, unsigned int const j, unsigned int const k );

/**
 * @brief Implements a cut cell criteria.
 * @param levelset The levelset field.
 * @param i The index in x-direction.
 * @param j The index in y-direction.
 * @param k The index in z-direction.
 * @return A bool indicating, whether the given index describes a cut cell.
 * @tparam Different cut cell criteria can b used. Template specifications implement the different criteria.
 */
template<CutCellCriteria>
bool IsCutCell( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()],
                unsigned int const i, unsigned int const j, unsigned int const k );

void ComputeInterfaceCurvature( Node const& node, double ( &curvature )[CC::TCX()][CC::TCY()][CC::TCZ()] );

void GetLevelsetAtSubcellCorners( double ( &subcell_corner_levelset )[subcell_box_size_x][subcell_box_size_y][subcell_box_size_z],
                                  double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                  unsigned int const i, unsigned int const j, unsigned int const k );

std::array<double, 3> GetNormal( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                 unsigned int const i, unsigned int const j, unsigned int const k, std::int8_t const material_sign = 1 );

#endif//GEOMETRY_CALCULATOR_H
