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
#ifndef LEVELSET_INITIALIZER_H
#define LEVELSET_INITIALIZER_H

#include "user_specifications/compile_time_constants.h"
#include "topology/node.h"

/**
 * @brief The LevelsetInitializer class defines an interface for different initialization methods for the levelset field. It contains the loop
 *        structure and the derived classes implement the levelset computation.
 * @note Each bounding box must fully enclose the negative field inside (e.g., if a plane is provided in 3D, the bounding box must extent to the full negative
 *       region and not only the interface).
 */
class LevelsetInitializer {
   // Member of this class
   std::vector<std::array<double, 6>> const bounding_boxes_;
   std::vector<MaterialName> const material_names_;
   double const dimensionalized_node_size_on_level_zero_;
   unsigned int const maximum_level_;

   // Functions for this class only
   bool IsPointInsideBoundingBox( std::array<double, 3> const& point ) const;

protected:
   // Protected member
   std::vector<std::string> const spatial_variable_names_ = { "x", "y", "z" };

   // Functions that need to be implemented by the derived classes
   virtual double ComputeSignedLevelsetValue( std::array<double, 3> const& point ) = 0;
   virtual std::string GetTypeLogData( unsigned int const indent ) const           = 0;

   // protected default constructor (can only be called from derived classes)
   explicit LevelsetInitializer( std::vector<std::array<double, 6>> const bounding_boxes,
                                 std::vector<MaterialName> const& material_names,
                                 double const node_size_on_level_zero,
                                 unsigned int const maximum_level );

public:
   virtual ~LevelsetInitializer()                    = default;
   LevelsetInitializer( LevelsetInitializer const& ) = delete;
   LevelsetInitializer& operator=( LevelsetInitializer const& ) = delete;
   LevelsetInitializer( LevelsetInitializer&& )                 = delete;
   LevelsetInitializer& operator=( LevelsetInitializer&& ) = delete;

   // Public function that can be called from outside
   void GetInitialLevelset( nid_t const node_id, double ( &initial_levelset )[CC::TCX()][CC::TCY()][CC::TCZ()] );
   std::vector<MaterialName> GetInitialMaterials( nid_t const node_id );
   std::string GetLogData( unsigned int const indent ) const;
};

#endif//LEVELSET_INITIALIZER_H
