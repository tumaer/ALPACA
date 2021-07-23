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

#include "initial_condition/levelset_initializer.h"

#include "utilities/mathematical_functions.h"
#include "utilities/buffer_operations.h"

/**
 * @brief Protected constructor that instantiates the base class of all levelset initializers.
 * @param bounding_boxes The bounding boxes inside which an interface can be found.
 * @param material_names Names of the materials.
 * @param node_size_on_level_zero Size of node on level zero.
 * @param maximum_level Maximum level of the simulation.
 */
LevelsetInitializer::LevelsetInitializer( std::vector<std::array<double, 6>> const bounding_boxes,
                                          std::vector<MaterialName> const& material_names,
                                          double const node_size_on_level_zero,
                                          unsigned int const maximum_level ) : bounding_boxes_( bounding_boxes ),
                                                                               material_names_( material_names ),
                                                                               dimensionalized_node_size_on_level_zero_( node_size_on_level_zero ),
                                                                               maximum_level_( maximum_level ) {
   /** Empty besides initializer list */
}

/**
 * @brief Determines whether a point is inside or outside the bounding boxes specified in the inputfile.
 * @param point Coordinates of the cell center.
 * @return True if point is inside, false otherwise.
 */
bool LevelsetInitializer::IsPointInsideBoundingBox( std::array<double, 3> const& point ) const {
   if( bounding_boxes_.size() == 0 ) {
      return true;
   }
   for( auto const& bounding_box : bounding_boxes_ ) {
      bool const x = point[0] > bounding_box[0] && point[0] < bounding_box[1];
      bool const y = CC::DIM() != Dimension::One ? point[1] > bounding_box[2] && point[1] < bounding_box[3] : true;
      bool const z = CC::DIM() == Dimension::Three ? point[2] > bounding_box[4] && point[2] < bounding_box[5] : true;
      if( x && y && z ) { return true; }
   }
   return false;
}

/**
 * @brief Gives the initial levelset value for a certain node ID.
 * @param node_id The node id for which the initial levelset should be obtained.
 * @param levelset The buffer, where the levelset values are written to.
 */
void LevelsetInitializer::GetInitialLevelset( nid_t const node_id, double ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {

   // For single materials always fill the levelset buffer with the cut off factor
   if( material_names_.size() == 1 ) {
      BO::SetSingleBuffer( levelset, CC::LSCOF() );
   } else {
      // get the origin of this node id
      std::array<double, 3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, dimensionalized_node_size_on_level_zero_ ) );
      double const cell_size             = CellSizeOfId( node_id, dimensionalized_node_size_on_level_zero_ );
      double const one_cell_size         = 1.0 / cell_size;

      // Coordinate vector that is filled on the fly
      std::array<double, 3> point = { 0.0, 0.0, 0.0 };

      // Loop through all cells to assign correct values to the buffer
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         point[0] = origin[0] + ( double( i ) - double( CC::FICX() ) + 0.5 ) * cell_size;
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            if constexpr( CC::DIM() != Dimension::One ) point[1] = origin[1] + ( double( j ) - double( CC::FICY() ) + 0.5 ) * cell_size;
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               if constexpr( CC::DIM() != Dimension::One ) point[2] = origin[2] + ( double( k ) - double( CC::FICZ() ) + 0.5 ) * cell_size;
               // Call the derived class function if the point is inside any bounding box
               levelset[i][j][k] = IsPointInsideBoundingBox( point ) ? ComputeSignedLevelsetValue( point ) * one_cell_size : CC::LSCOF();
            }
         }
      }
   }
}

/**
 * @brief Gives the initial materials for a certain node ID.
 * @param node_id The node id for which the initial materials should be obtained.
 * @return The materials contained in this node.
 */
std::vector<MaterialName> LevelsetInitializer::GetInitialMaterials( nid_t const node_id ) {

   // For single materials always return the material immediately
   if( material_names_.size() == 1 ) {
      return material_names_;
   }

   // get the origin of this node id
   std::array<double, 3> const origin      = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, dimensionalized_node_size_on_level_zero_ ) );
   double const cell_size_on_maximum_level = CellSizeOfLevel( dimensionalized_node_size_on_level_zero_, maximum_level_ );

   // Level factors for difference between current node level and maximum level of simulation
   // bit shift is of type "( unsigned? ) int"
   unsigned int const level_factor_x = ( 1 << ( maximum_level_ - LevelOfNode( node_id ) ) );
   unsigned int const level_factor_y = CC::DIM() != Dimension::One ? level_factor_x : 1;
   unsigned int const level_factor_z = CC::DIM() == Dimension::Three ? level_factor_x : 1;

   // Vector that is returned (yet not known which size it has)
   std::vector<MaterialName> materials_contained;
   // Vector containing flags for each material if it was already added for this node
   std::vector<bool> materials_present( material_names_.size(), false );

   // // Temporary point that is filled
   std::array<double, 3> point = { 0.0, 0.0, 0.0 };

   // In the following loop we run over all cells on the finest level that would cover the same region as the node of interest including its halo region
   for( unsigned int i = 0; i < level_factor_x * CC::TCX(); ++i ) {
      point[0] = origin[0] + ( double( i ) - double( level_factor_x * CC::FICX() ) + 0.5 ) * cell_size_on_maximum_level;
      for( unsigned int j = 0; j < level_factor_y * CC::TCY(); ++j ) {
         if constexpr( CC::DIM() != Dimension::One ) point[1] = origin[1] + ( double( j ) - double( level_factor_y * CC::FICY() ) + 0.5 ) * cell_size_on_maximum_level;
         for( unsigned int k = 0; k < level_factor_z * CC::TCZ(); ++k ) {
            if constexpr( CC::DIM() == Dimension::Three ) point[2] = origin[2] + ( double( k ) - double( level_factor_z * CC::FICZ() ) + 0.5 ) * cell_size_on_maximum_level;
            // Call the derived class function if the point is inside any bounding box
            int const levelset_sign = IsPointInsideBoundingBox( point ) ? Signum( ComputeSignedLevelsetValue( point ) ) : 1;
            // Add the material depending on the sign. For zero immediately return both materials
            if( levelset_sign < 0 && !materials_present[1] ) {
               materials_present[1] = true;
               materials_contained.push_back( material_names_[1] );
            } else if( levelset_sign > 0 && !materials_present[0] ) {
               materials_present[0] = true;
               materials_contained.push_back( material_names_[0] );
            }
            // If all materials are contained, we can already return
            if( material_names_.size() == materials_contained.size() ) return materials_contained;
         }
      }
   }

   // Return the contained material
   return materials_contained;
}

/**
 * @brief Gives the data for logging for this appropriate class.
 * @param indent Number of white spaces used at the beginning of each line for the logging information.
 * @return string with logging information.
 */
std::string LevelsetInitializer::GetLogData( unsigned int const indent ) const {
   // The log string that is created (call the derived class function )
   std::string log_string( GetTypeLogData( indent ) );
   log_string += "\nBounding boxes:\n";
   // Loop through all bounding boxes and add the successively
   if( bounding_boxes_.size() == 0 ) {
      log_string += "No boxes provided";
   } else {
      for( unsigned int index = 0; index < bounding_boxes_.size(); index++ ) {
         // Convert all bounding box values to scientific notation
         std::array<std::string, 6> bb_string;
         std::transform( std::cbegin( bounding_boxes_[index] ), std::cend( bounding_boxes_[index] ), std::begin( bb_string ),
                         []( double const lim ) { return StringOperations::ToScientificNotationString( lim, 2 ); } );
         // Create the log string
         log_string += StringOperations::Indent( indent + 2 ) + "Box " + std::to_string( index + 1 ) + "\n";
         for( unsigned int dim = 0; dim < DTI( CC::DIM() ); ++dim ) {
            log_string += StringOperations::Indent( indent + 4 ) + spatial_variable_names_[dim] + ": " + bb_string[2 * dim] + " " + bb_string[2 * dim + 1] + "\n";
         }
      }
   }
   return log_string;
}
