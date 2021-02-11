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
#include "parametric_levelset_initializer.h"

#include <stdexcept>
#include <cmath>

#include "utilities/mathematical_functions.h"
#include "utilities/helper_functions.h"
#include "utilities/string_operations.h"
#include "utilities/vector_utilities.h"
#include "user_expression.h"

/**
 * @brief Computes all interface coordinates for a given expression.
 * @param parametric_expression The expression specifying the interface coordinates.
 * @param parametric_variables The parametric variables to be used to evaluate the expression.
 * @return All interface coordinates.
 */
std::vector<std::array<double, 3>> ComputeInterfaceCoordinates( std::string const& parametric_expression,
                                                                std::vector<std::string> const& spatial_variable_names,
                                                                std::array<ParametricVariable, 2> const& parametric_variables ) {

   // Extract values for better handling
   std::vector<std::string> const variables_names      = { parametric_variables[0].name, parametric_variables[1].name };
   std::array<double, 2> const start_values            = { parametric_variables[0].start, parametric_variables[1].start };
   std::array<double, 2> const delta_increments        = { parametric_variables[0].delta, parametric_variables[1].delta };
   std::array<std::uint64_t, 2> const number_of_points = { parametric_variables[0].points, parametric_variables[1].points };
   std::uint64_t const total_number_of_points          = number_of_points[0] * number_of_points[1];

   // Define the total number of points and the interface coordinates vector
   std::vector<std::array<double, 3>> interface_coordinates;
   interface_coordinates.reserve( total_number_of_points );

   // Create the parametric expression with two running coordinates
   std::vector<double> parameteric_point( 2, 0.0 );
   UserExpression const parameteric_expression( UserExpression( parametric_expression, spatial_variable_names, variables_names, parameteric_point ) );

   std::array<double, 3> interface_point = { 0.0, 0.0, 0.0 };

   for( std::uint64_t i = 0; i < number_of_points[0]; ++i ) {
      parameteric_point[0] = start_values[0] + double( i ) * delta_increments[0];
      for( std::uint64_t j = 0; j < number_of_points[1]; ++j ) {
         parameteric_point[1] = start_values[1] + double( j ) * delta_increments[1];
         std::transform( std::cbegin( spatial_variable_names ), std::cend( spatial_variable_names ), std::begin( interface_point ),
                         [&parameteric_expression]( std::string const variable ) { return parameteric_expression.GetValue( variable ); } );
         interface_coordinates.push_back( interface_point );
      }
   }
   return interface_coordinates;
}

/**
 * @brief Constructs a parametric levelset initializer with levelset initialization parameters given as input.
 * @param parameteric_expression_string The expression string for the parametric interface.
 * @param parametric_variables The parameteric variables for the interface computation.
 * @param ref_point_of_negative_levelset The reference point in the negative material region.
 * @param bounding_boxes The bounding boxes inside which an interface can be found.
 * @param material_names Names of the materials.
 * @param node_size_on_level_zero Size of node on level zero.
 * @param maximum_level Maximum level of the simulation.
 */
ParametricLevelsetInitializer::ParametricLevelsetInitializer( std::string const& parameteric_expression_string,
                                                              std::array<ParametricVariable, 2> const& parametric_variables,
                                                              std::array<double, 3> const& ref_point_of_negative_levelset,
                                                              std::vector<std::array<double, 6>> const& bounding_boxes,
                                                              std::vector<MaterialName> const& material_names,
                                                              double const node_size_on_level_zero,
                                                              unsigned int const maximum_level ) : LevelsetInitializer( bounding_boxes, material_names, node_size_on_level_zero, maximum_level ),
                                                                                                   parametric_levelset_expession_( parameteric_expression_string ),
                                                                                                   parametric_variables_( parametric_variables ),
                                                                                                   ref_point_of_negative_levelset_( ref_point_of_negative_levelset ),
                                                                                                   interface_coordinates_( ComputeInterfaceCoordinates( parametric_levelset_expession_, spatial_variable_names_, parametric_variables_ ) ) {
   /* Empty besides initializer list*/
}

/**
 * @brief Finds the minimum distance and corresponding point on the interface.
 * @param point The point for which the minimum distance should be obtained.
 * @return The point and distance.
 */
std::pair<std::array<double, 3>, double> ParametricLevelsetInitializer::FindInterfacePointWithMinimumDistance( std::array<double, 3> const& point ) const {
   // Compute the distances between the interface points and the current point
   std::vector<double> distances( interface_coordinates_.size(), std::numeric_limits<double>::max() );
   std::transform( std::cbegin( interface_coordinates_ ), std::cend( interface_coordinates_ ), std::begin( distances ),
                   [&point]( std::array<double, 3> const& interface_point ) { return VU::Distance( point, interface_point ); } );
   // Get the minimum distance iterator
   auto const nearest_point_it = std::min_element( std::cbegin( distances ), std::cend( distances ) );
   // Return the distance and the point
   return std::make_pair( interface_coordinates_[nearest_point_it - std::cbegin( distances )], *nearest_point_it );
}

/**
 * @brief Gives the levelset sign for a point compared to its nearest interface point.
 * @param point The point for which the sign should be obtained.
 * @param interface_point The point of the interface with the nearest distance to point.
 * @return The levelset sign.
 */
int ParametricLevelsetInitializer::GetLevelsetSign( std::array<double, 3> const& point, std::array<double, 3> const& interface_point ) const {
   // Compute the differences between the point and the center and the interface and the center
   double const interface_to_ref_distance = VU::Distance( interface_point, ref_point_of_negative_levelset_ );
   double const point_to_ref_distance     = VU::Distance( point, ref_point_of_negative_levelset_ );
   // Return the distance with the appropriate sign
   return point_to_ref_distance < interface_to_ref_distance ? -1 : 1;
}

/**
 * @brief Computes the level value and appropriate sign for the given point.
 * @param point The point for which the levelset value and sign should be obtained.
 * @return The signed levelset value.
 */
double ParametricLevelsetInitializer::ComputeSignedLevelsetValue( std::array<double, 3> const& point ) {
   std::pair<std::array<double, 3>, double> const interface_point_and_distance = FindInterfacePointWithMinimumDistance( point );
   return interface_point_and_distance.second * static_cast<double>( GetLevelsetSign( point, interface_point_and_distance.first ) );
}

/**
 * @brief Gives the data for logging for this appropriate class.
 * @param indent Number of white spaces used at the beginning of each line for the logging information.
 * @return string with logging information.
 */
std::string ParametricLevelsetInitializer::GetTypeLogData( unsigned int const indent ) const {
   // string that is returned
   std::string log_string;
   // Name of the levelset initializer
   log_string += StringOperations::Indent( indent ) + "Type     : Parametric\n";
   // Function
   log_string += StringOperations::Indent( indent ) + "Function :\n";
   log_string += StringOperations::Indent( indent + 2 ) + parametric_levelset_expession_ + "\n\n";
   // Reference point
   log_string += "Reference point negative levelset:\n";
   for( unsigned int dim = 0; dim < DTI( CC::DIM() ); ++dim ) {
      log_string += StringOperations::Indent( indent + 2 ) + spatial_variable_names_[dim] + ": " +
                    StringOperations::ToScientificNotationString( ref_point_of_negative_levelset_[dim], 4 ) + "\n";
   }
   // Parameters
   log_string += "Parameters:\n";
   for( auto const& parameter : parametric_variables_ ) {
      if( !parameter.name.empty() ) {
         log_string += parameter.GetLogData( indent ) + "\n";
      }
   }

   return log_string;
}
