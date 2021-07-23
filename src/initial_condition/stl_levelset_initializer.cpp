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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#include "stl_levelset_initializer.h"

/**
 * @brief Constructs a stl levelset initializer with levelset initialization parameters given as input.
 * @param dimensional_levelset_initializer_data Map containing all data for the levelset initializer.
 * @param bounding_boxes The bounding boxes inside which an interface can be found.
 * @param material_names Names of the materials.
 * @param node_size_on_level_zero Size of node on level zero.
 * @param maximum_level Maximum level of the simulation.
 */
StlLevelsetInitializer::StlLevelsetInitializer( std::string const& stl_filename,
                                                std::vector<std::array<double, 6>> const& bounding_boxes,
                                                std::vector<MaterialName> const& material_names,
                                                double const node_size_on_level_zero,
                                                unsigned int const maximum_level ) : LevelsetInitializer( bounding_boxes, material_names, node_size_on_level_zero, maximum_level ),
                                                                                     stl_filename_( stl_filename ),
                                                                                     stl_triangles_( StlUtilities::ReadStl( stl_filename_ ) ) {
   /* Empty besides initializer list*/
}

/**
 * @brief Computes the level value and appropriate sign for the given point.
 * @param point The point for which the levelset value and sign should be obtained.
 * @return The signed levelset value.
 */
double StlLevelsetInitializer::ComputeSignedLevelsetValue( std::array<double, 3> const& point ) {
   double levelset = std::numeric_limits<double>::max();
   for( StlUtilities::Triangle const& triangle : stl_triangles_ ) {
      StlUtilities::Voxelization( triangle, point, levelset );
   }
   return levelset;
}

/**
 * @brief Gives the data for logging for this appropriate class.
 * @param indent Number of white spaces used at the beginning of each line for the logging information.
 * @return string with logging information.
 */
std::string StlLevelsetInitializer::GetTypeLogData( unsigned int const indent ) const {
   // string that is returned
   std::string log_string;
   // Name of the levelset initializer
   log_string += StringOperations::Indent( indent ) + "Type               : STL\n";
   // Function
   log_string += StringOperations::Indent( indent ) + "Filename           : " + stl_filename_ + "\n";
   // Number of triangles
   log_string += StringOperations::Indent( indent ) + "Number of triangles: " + std::to_string( stl_triangles_.size() ) + "\n";

   return log_string;
}
