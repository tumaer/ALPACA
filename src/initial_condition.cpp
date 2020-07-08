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
#include "initial_condition.h"

#include <stdexcept>
#include "topology/id_information.h"

/**
 * @brief Default constructor for class Converting the user-expressions for the initial conditions.
 * @param material_initial_expressions Initial expressions for all materials contained in the simulation.
 * @param levelset_initial_expressions Initial expressions for all levelset fields contained in the simulation
 *        (All must be referenced as a function between material-N and material-1).
 * @param material_names List of all material names contained in the simulation.
 * @param variable_names_prime_states List of all names for the prime states that are read for the initial conditions.
 * @param variable_name_levelset Name of the variable used to specify levelset field.
 * @param dimensionalized_node_size_on_level_zero Size of a node on level zero (dimensionalized form).
 * @param maximum_level Maximum level contained in the simulation.
 * @param unit_handler Instance to provide (non-)dimensionalization operations.
 */
InitialCondition::InitialCondition( std::vector<std::string> const& material_initial_expressions,
                                    std::vector<std::string> const& levelset_initial_expressions,
                                    std::vector<MaterialName> const& material_names,
                                    std::vector<std::string> const& variable_names_prime_states,
                                    std::string const& variable_name_levelset,
                                    double const dimensionalized_node_size_on_level_zero,
                                    unsigned int const maximum_level,
                                    UnitHandler const& unit_handler ) :// Start initializer list
                                                                        unit_handler_( unit_handler ),
                                                                        material_initial_expressions_( material_initial_expressions ),
                                                                        levelset_initial_expressions_( levelset_initial_expressions ),
                                                                        material_names_( material_names ),
                                                                        variable_names_prime_states_( variable_names_prime_states ),
                                                                        variable_name_levelset_( variable_name_levelset ),
                                                                        dimensionalized_node_size_on_level_zero_( dimensionalized_node_size_on_level_zero ),
                                                                        maximum_level_( maximum_level ) {
   /** Empty besides initializer list */
}

/**
 * @brief Compiles the given expression such that it is capable to give the desired material variables.
 * @param expression The expression in original text form.
 * @param x,y,z Reference to spatial variable.
 * @return The compiled expression ready to be evaluated.
 *
 * @note The return type must be a pointer to allow the creation of UserExpressions inside of a vector (required for the initial materials). This is due to
 *       the deleted move operators of the user_expression library
 */
std::unique_ptr<UserExpression const> InitialCondition::CreateInputExpression( std::string const& expression,
                                                                               std::vector<std::string> const& variables_out,
                                                                               double& x, double& y, double& z ) const {
   std::vector<std::tuple<std::string, double&>> variables_in;
   variables_in.push_back( std::make_tuple( std::string( variable_name_x_ ), std::ref( x ) ) );
   variables_in.push_back( std::make_tuple( std::string( variable_name_y_ ), std::ref( y ) ) );
   variables_in.push_back( std::make_tuple( std::string( variable_name_z_ ), std::ref( z ) ) );

   return std::make_unique<UserExpression const>( expression, variables_in, variables_out );
}

/**
 * @brief Gives the initial density at the provided location for the given material.
 * @param node_id The id of the node to be initialized.
 * @param material The material in the block to be filled with the returned data.
 * @param initial_values Reference to array holding the resulting density. Indirect return value.
 */
void InitialCondition::GetInitialPrimeStates( std::uint64_t const node_id,
                                              MaterialName const material,
                                              double ( &initial_values )[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()] ) const {

   // get the origin of this node id
   std::array<double, 3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, dimensionalized_node_size_on_level_zero_ ) );
   // cell_size on level
   double const cell_size = dimensionalized_node_size_on_level_zero_ / double( CC::ICX() ) / double( 1 << LevelOfNode( node_id ) );

   // Obtain the correct initial condition string for the material
   std::string const initial_condition_input( material_initial_expressions_[MTI( material )] );

   // TP here we need non-const variables as UserExpression stores references to them in order to reflect changes
   double running_x;
   double running_y;
   double running_z;

   // create the expression
   std::unique_ptr<UserExpression const> input_expression( CreateInputExpression( initial_condition_input, variable_names_prime_states_, running_x, running_y, running_z ) );
   // Loop through all cells to assign correct values to the buffer
   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      running_x = origin[0] + ( double( i ) + 0.5 ) * cell_size;
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         running_y = CC::DIM() != Dimension::One ? origin[1] + ( double( j ) + 0.5 ) * cell_size : 0.0;
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            running_z = CC::DIM() == Dimension::Three ? origin[2] + ( double( k ) + 0.5 ) * cell_size : 0.0;
            for( PrimeState const p : MF::ASOP() ) {
               // If the variable name is not empty obtain value from expression.
               if( !variable_names_prime_states_[PTI( p )].empty() ) {
                  initial_values[PTI( p )][i][j][k] = unit_handler_.NonDimensionalizeValue( input_expression->GetValue( variable_names_prime_states_[PTI( p )] ), MF::FieldUnit( p ) );
               } else {// Otherwise set zero value
                  initial_values[PTI( p )][i][j][k] = 0.0;
               }
            }
         }
      }
   }
}

/**
 * @brief Gives the initial levelset for the node with the given id.
 * @param node_id The id of the node to be initialized.
 * @param initial_levelset Indirect return parameter for the determined levelset.
 */
void InitialCondition::GetInitialLevelset( std::uint64_t const node_id, double ( &initial_levelset )[CC::TCX()][CC::TCY()][CC::TCZ()] ) const {
   // get the origin of this node id
   std::array<double, 3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, dimensionalized_node_size_on_level_zero_ ) );
   // cell_size on level zero divided by 2^level
   double const cell_size     = dimensionalized_node_size_on_level_zero_ / double( CC::ICX() ) / double( 1 << LevelOfNode( node_id ) );
   double const one_cell_size = 1.0 / cell_size;

   // NOTE: This must be changed for Multi-material approach (all initial levelset should be referenced to material1 as done in GetInitialMaterials)
   std::string initial_levelset_input = levelset_initial_expressions_[0];

   // TP here we need non-const variables as UserExpression stores references to them in order to reflect changes
   double running_x;
   double running_y;
   double running_z;

   // create the expression
   std::unique_ptr<UserExpression const> levelset_expression( CreateInputExpression( initial_levelset_input, { variable_name_levelset_ }, running_x, running_y, running_z ) );
   // Loop through all cells to assign correct values to tje buffer
   for( unsigned int i = 0; i < CC::TCX(); ++i ) {
      running_x = origin[0] + ( double( i ) - double( CC::FICX() ) + 0.5 ) * cell_size;
      for( unsigned int j = 0; j < CC::TCY(); ++j ) {
         running_y = CC::DIM() != Dimension::One ? origin[1] + ( double( j ) - double( CC::FICY() ) + 0.5 ) * cell_size : 0.0;
         for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
            running_z = CC::DIM() == Dimension::Three ? origin[2] + ( double( k ) - double( CC::FICZ() ) + 0.5 ) * cell_size : 0.0;
            // Non-dimensionalize value with cell size
            initial_levelset[i][j][k] = levelset_expression->GetValue( variable_name_levelset_ ) * one_cell_size;
         }
      }
   }
}

/**
 * @brief Gives the initial materials for the node with the given id.
 * @param node_id The id of the node to be initialized.
 * @return A list of the materials present in the considered node.
 */
std::vector<MaterialName> InitialCondition::GetInitialMaterials( std::uint64_t const node_id ) const {
   // get the origin of this node id
   std::array<double, 3> const origin = DomainCoordinatesOfId( node_id, DomainSizeOfId( node_id, dimensionalized_node_size_on_level_zero_ ) );
   // bit shift is of type "( unsigned? ) int"
   unsigned int const level_factor = ( 1 << ( maximum_level_ - LevelOfNode( node_id ) ) );
   // cell size on maximum level
   double const cell_size_on_maximum_level = dimensionalized_node_size_on_level_zero_ / double( CC::ICX() ) / double( 1 << maximum_level_ );

   //************************************************************************************************************
   // Vector that is returned (yet not known which size it has)
   std::vector<MaterialName> materials_contained;
   std::vector<bool> material_is_contained( material_names_.size() );

   // Single material can be returned immediately
   if( material_names_.size() == 1 ) {
      materials_contained.push_back( material_names_[0] );
      return materials_contained;
   } else {

      // TP here we need non-const variables as UserExpression stores references to them in order to reflect changes
      double running_x;
      double running_y;
      double running_z;

      // create all levelset expressions (total number is number of materials - 1). All materials are referred to material One which has a positive value.
      // All other take a negative value
      std::vector<std::unique_ptr<UserExpression const>> levelset_expressions;
      levelset_expressions.reserve( levelset_initial_expressions_.size() );
      for( auto const& input_expression : levelset_initial_expressions_ ) {
         levelset_expressions.emplace_back( CreateInputExpression( input_expression, { variable_name_levelset_ }, running_x, running_y, running_z ) );
      }

      // Level factors for difference between current node level and maximum level of simulation
      unsigned int const level_factor_y = CC::DIM() != Dimension::One ? level_factor : 1;
      unsigned int const level_factor_z = CC::DIM() == Dimension::Three ? level_factor : 1;

      // variables required to identify correct assignment of present materials
      double levelset_value;
      // Vector containing flags for each material if it was already added for this node
      std::vector<bool> materials_present( material_names_.size(), false );

      // In the following loop we run over all cells on the finest level that would cover the same region as the node of interest including its halo region
      for( unsigned int i = 0; i < level_factor * CC::TCX(); ++i ) {
         running_x = origin[0] + ( double( i ) - double( level_factor * CC::FICX() ) + 0.5 ) * cell_size_on_maximum_level;
         for( unsigned int j = 0; j < level_factor_y * CC::TCY(); ++j ) {
            running_y = CC::DIM() != Dimension::One ? origin[1] + ( double( j ) - double( level_factor_y * CC::FICY() ) + 0.5 ) * cell_size_on_maximum_level : 0.0;
            for( unsigned int k = 0; k < level_factor_z * CC::TCZ(); ++k ) {
               running_z = CC::DIM() == Dimension::Three ? origin[2] + ( double( k ) - double( level_factor_z * CC::FICZ() ) + 0.5 ) * cell_size_on_maximum_level : 0.0;

               // Flag that a cell already contains another negative material
               bool cell_contains_negative_material = false;
               // counter for the levelset expressions (start at 1 since MaterialOne does not have an expression)
               unsigned int levelset_counter = 1;

               // Get all levelset values at the current cell
               for( auto const& levelset_expression : levelset_expressions ) {
                  levelset_value = levelset_expression->GetValue( variable_name_levelset_ );
                  if( levelset_value < 0.0 ) {
                     // If the flag, that a negative material was already found in this cell, is true
                     // throw error since two negative materials at the same time are not allowed
                     if( cell_contains_negative_material ) {
                        throw std::logic_error( "Error! At least two materials have negative levelset values in the same cell. "
                                                "Check initial levelset conditions for all materials!" );
                     } else {// Otherwise take the negative and assign it
                        cell_contains_negative_material = true;
                        // Only set the flag and append to final vector if not already done for this material
                        if( !materials_present[levelset_counter] ) {
                           materials_present[levelset_counter] = true;
                           materials_contained.push_back( material_names_[levelset_counter] );
                        }
                     }
                  }
                  // increment counter for next expression
                  levelset_counter++;
               }

               // If no cell has a negative value append the positive (MaterialOne) if not already done
               if( !cell_contains_negative_material && !materials_present[0] ) {
                  materials_present[0] = true;
                  materials_contained.push_back( material_names_[0] );
               }

               // If all materials are contained, we can already return
               if( material_names_.size() == materials_contained.size() ) return materials_contained;
            }
         }
      }

      return materials_contained;
   }
}