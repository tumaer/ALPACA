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
#include "prime_state_initializer.h"

#include "topology/id_information.h"
#include "user_expression.h"

/**
 * @brief Constructs the prime state initializer to evaluate prime state values from input expressions.
 * @param prime_state_expression_strings The expression string for all materials contained in the simulation.
 * @param prime_state_variable_names List of all names for the prime states that are read for the initial conditions.
 * @param dimensionalized_node_size_on_level_zero Size of a node on level zero (dimensionalized form).
 * @param unit_handler Instance to provide (non-)dimensionalization operations.
 */
PrimeStateInitializer::PrimeStateInitializer( std::vector<std::string> const& prime_state_expression_strings,
                                              std::vector<std::string> const& prime_state_variable_names,
                                              double const dimensionalized_node_size_on_level_zero,
                                              UnitHandler const& unit_handler ) : unit_handler_( unit_handler ),
                                                                                  prime_state_expression_strings_( prime_state_expression_strings ),
                                                                                  prime_state_variable_names_( prime_state_variable_names ),
                                                                                  dimensionalized_node_size_on_level_zero_( dimensionalized_node_size_on_level_zero ) {
   /** Empty besides initializer list */
}

/**
 * @brief Gives the initial prime states for the given material.
 * @param node_id The id of the node to be initialized.
 * @param material The material in the block to be filled with the returned data.
 * @param prime_state_buffer Reference to array holding the resulting prime states. Indirect return value.
 */
void PrimeStateInitializer::GetInitialPrimeStates( nid_t const node_id,
                                                   MaterialName const material,
                                                   double ( &prime_state_buffer )[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()] ) const {

   // get the origin of this node id
   std::array<double, 3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, dimensionalized_node_size_on_level_zero_ ) );
   double const cell_size             = CellSizeOfId( node_id, dimensionalized_node_size_on_level_zero_ );

   // create the expression
   std::vector<double> cell_center_point( 3, 0.0 );
   UserExpression const prime_state_expression( UserExpression( prime_state_expression_strings_[MTI( material )], prime_state_variable_names_,
                                                                spatial_variable_names_, cell_center_point ) );

   // Loop through all cells to assign correct values to the buffer
   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      cell_center_point[0] = origin[0] + ( double( i ) + 0.5 ) * cell_size;
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         if constexpr( CC::DIM() != Dimension::One ) cell_center_point[1] = origin[1] + ( double( j ) + 0.5 ) * cell_size;
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            if constexpr( CC::DIM() == Dimension::Three ) cell_center_point[2] = origin[2] + ( double( k ) + 0.5 ) * cell_size;
            for( PrimeState const p : MF::ASOP() ) {
               // If the variable name is not empty obtain value from expression.
               if( !prime_state_variable_names_[PTI( p )].empty() ) {
                  prime_state_buffer[PTI( p )][i][j][k] = unit_handler_.NonDimensionalizeValue( prime_state_expression.GetValue( prime_state_variable_names_[PTI( p )] ), MF::FieldUnit( p ) );
               } else {// Otherwise set zero value
                  prime_state_buffer[PTI( p )][i][j][k] = 0.0;
               }
            }
         }
      }
   }
}