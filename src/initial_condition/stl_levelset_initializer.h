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
#ifndef STL_LEVELSET_INITIALIZER_H
#define STL_LEVELSET_INITIALIZER_H

#include "initial_condition/levelset_initializer.h"
#include "input_output/utilities/stl_utilities.h"

/**
 * @brief The StlLevelsetInitializer class allows for a levelset initialization based on a provided STL file. It contains the levelset computation,
 *        the loop structure is in the base class.
  * @note The number and size of the triangles, the surface is composed of, must be appropriate for the used mesh size. In case smaller cell sizes are used,
  *       the number of triangles must be increased and the size reduced. The inversely holds for coarser meshes.
 */
class StlLevelsetInitializer : public LevelsetInitializer {
   // Member variables of this class only
   std::string const stl_filename_;

   std::vector<StlUtilities::Triangle> stl_triangles_;

   // Functions required from base class
   double ComputeSignedLevelsetValue( std::array<double, 3> const& point ) override;
   std::string GetTypeLogData( unsigned int const indent ) const override;

public:
   StlLevelsetInitializer() = delete;
   explicit StlLevelsetInitializer( std::string const& stl_filename,
                                    std::vector<std::array<double, 6>> const& bounding_boxes,
                                    std::vector<MaterialName> const& material_names,
                                    double const node_size_on_level_zero,
                                    unsigned int const maximum_level );
   virtual ~StlLevelsetInitializer()                       = default;
   StlLevelsetInitializer( StlLevelsetInitializer const& ) = delete;
   StlLevelsetInitializer& operator=( StlLevelsetInitializer const& ) = delete;
   StlLevelsetInitializer( StlLevelsetInitializer&& )                 = delete;
   StlLevelsetInitializer& operator=( StlLevelsetInitializer&& ) = delete;
};

#endif//STL_LEVELSET_INITIALIZER_H
