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
#ifndef PARAMETRIC_LEVELSET_INITIALIZER_H
#define PARAMETRIC_LEVELSET_INITIALIZER_H

#include "initial_condition/levelset_initializer.h"
#include "initial_condition/parametric_variable.h"

/**
 * @brief The ParametricLevelsetInitializer class allows for a levelset initialization based on parametric functions. It contains the levelset computation,
 *        the loop structure is in the base class.
 * @note Uses the C++ Mathematical Expression Toolkit Library by Arash Partow, see respective files for License and Copyright information.
 * @note For the parametric levelset some things need to be considered:
 *       1. The number of points for each parametric variable must be large enough to allow proper sampling compared to the cell size.
 *       2. The reference point can be anywhere in the negative material region (for open interfaces, such as quasi 1D in 3D, the reference point should be
 *          placed outside of the domain to ensure that no cells are behind the reference point).
 *       3. The parametric function must be also valid in a small region outside of the domain (halo cells).
 *       4. The parametric function does not contain any undercuts.
 */
class ParametricLevelsetInitializer : public LevelsetInitializer {

   // Data obtained through the constructor
   std::string const parametric_levelset_expession_;
   std::array<ParametricVariable, 2> const parametric_variables_;
   std::array<double, 3> const ref_point_of_negative_levelset_;
   std::vector<std::array<double, 3>> interface_coordinates_;

   // Functions to compute the levelset values
   std::pair<std::array<double, 3>, double> FindInterfacePointWithMinimumDistance( std::array<double, 3> const& point ) const;
   int GetLevelsetSign( std::array<double, 3> const& point, std::array<double, 3> const& interface_point ) const;

   // Functions required from base class
   double ComputeSignedLevelsetValue( std::array<double, 3> const& point ) override;
   std::string GetTypeLogData( unsigned int const indent ) const override;

public:
   ParametricLevelsetInitializer() = delete;
   explicit ParametricLevelsetInitializer( std::string const& parameteric_expression_string,
                                           std::array<ParametricVariable, 2> const& parametric_variables,
                                           std::array<double, 3> const& ref_point_of_positive_levelset,
                                           std::vector<std::array<double, 6>> const& bounding_boxes,
                                           std::vector<MaterialName> const& material_names,
                                           double const node_size_on_level_zero,
                                           unsigned int const maximum_level );
   virtual ~ParametricLevelsetInitializer()                              = default;
   ParametricLevelsetInitializer( ParametricLevelsetInitializer const& ) = delete;
   ParametricLevelsetInitializer& operator=( ParametricLevelsetInitializer const& ) = delete;
   ParametricLevelsetInitializer( ParametricLevelsetInitializer&& )                 = delete;
   ParametricLevelsetInitializer& operator=( ParametricLevelsetInitializer&& ) = delete;
};

// Factory functions (outside for testing)
std::vector<std::array<double, 3>> ComputeInterfaceCoordinates( std::string const& parametric_expression,
                                                                std::vector<std::string> const& spatial_variable_names,
                                                                std::array<ParametricVariable, 2> const& parametric_variables );

#endif//PARAMETRIC_LEVELSET_INITIALIZER_H
