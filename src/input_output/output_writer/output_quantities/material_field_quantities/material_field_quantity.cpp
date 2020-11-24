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
#include "input_output/output_writer/output_quantities/material_field_quantities/material_field_quantity.h"
#include "input_output/output_writer/xdmf_utilities.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief constructor to create material field output.
 * @param unit_handler Unit handler class for dimensionalization.
 * @param material_manager Material manager for accessing material data.
 * @param quantity_name The name of the material field quantity used for the output.
 * @param output_flags Flags for which output type an output is written (0: standard, 1: interface, 2:debug).
 * @param quantity_data Data struct that contains all relevant information of the material quantity.
 * @param buffer_type Conservative buffert type identifier (default: Average).
 *
 * @note The quantity_name is written in the file (seen in ParaView cell_data list).
 */
MaterialFieldQuantity::MaterialFieldQuantity( UnitHandler const& unit_handler,
                                              MaterialManager const& material_manager,
                                              std::string const& quantity_name,
                                              std::array<bool, 3> const output_flags,
                                              MaterialFieldQuantityData const& quantity_data,
                                              ConservativeBufferType const buffer_type ) :// Start initializer list
                                                                                           OutputQuantity( unit_handler, material_manager, quantity_name + SuffixOfConservativeBufferType( buffer_type ), output_flags, quantity_data.GetDimensions() ),
                                                                                           quantity_data_( quantity_data ),
                                                                                           buffer_type_( buffer_type ) {
   /** Empty besides initializer list and base class constructor call */
}

/**
 * @brief See base class definition.
 */
void MaterialFieldQuantity::DoComputeCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int& cell_data_counter ) const {

   // extract the correct field type
   MaterialFieldType const field_type = quantity_data_.field_type_;

   // Change calls for levelset and non-levelset nodes
   if( node.HasLevelset() ) {

      std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
      double const( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()]            = node.GetInterfaceBlock().GetBaseBuffer( InterfaceDescription::Levelset );

      // local counter
      unsigned long long int local_counter = 0;

      //Loop through all components
      for( std::size_t component = 0; component < quantity_data_.field_indices_.size(); component++ ) {
         // Get the correct field index
         unsigned int const field_index = quantity_data_.field_indices_[component];
         // Get buffers of both materials and other specifications
         double const( &positive_field_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] =
               node.GetPhaseByMaterial( MaterialSignCapsule::PositiveMaterial() ).GetFieldBuffer( field_type, field_index, buffer_type_ );
         double const( &negative_field_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] =
               node.GetPhaseByMaterial( MaterialSignCapsule::NegativeMaterial() ).GetFieldBuffer( field_type, field_index, buffer_type_ );
         // Dimensionalization factor for re-dimensionalization of variables
         double const dimensionalization_factor = unit_handler_.DimensionalizeValue( 1.0, MF::FieldUnit( field_type, field_index ) );

         // set the local counter on original cell_data_counter + the component
         local_counter = cell_data_counter + component;

         // Loop through all internal cells of block
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
                  // Use negative material if interface tags are negative or in new cut cell band and negative levelset
                  if( interface_tags[i][j][k] < 0 || ( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell ) && levelset[i][j][k] < 0.0 ) ) {
                     cell_data[local_counter] = negative_field_buffer[i][j][k] * dimensionalization_factor;
                  }
                  // otherwise positive
                  else {
                     cell_data[local_counter] = positive_field_buffer[i][j][k] * dimensionalization_factor;
                  }

                  local_counter += quantity_data_.field_indices_.size();
               }
            }
         }
      }

   } else {
      // Obtain material and corresponding field buffer
      // No interface node -> interface tags/material is the same everywhere
      MaterialName const material = node.GetSinglePhaseMaterial();

      // local counter
      unsigned long long int local_counter = 0;

      //Loop through all components
      for( std::size_t component = 0; component < quantity_data_.field_indices_.size(); component++ ) {
         // Get the correct buffer and dimensionalization factor
         unsigned int const field_index                                 = quantity_data_.field_indices_[component];
         double const( &field_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( material ).GetFieldBuffer( field_type, field_index, buffer_type_ );
         double const dimensionalization_factor                         = unit_handler_.DimensionalizeValue( 1.0, MF::FieldUnit( field_type, field_index ) );

         // set the local counter on original cell_data_counter + the component
         local_counter = cell_data_counter + component;

         // Loop through all internal cells in the block
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
                  // add the dimensionalized value
                  cell_data[local_counter] = field_buffer[i][j][k] * dimensionalization_factor;

                  local_counter += quantity_data_.field_indices_.size();
               }
            }
         }
      }
   }
}

/**
 * @brief See base class definition.
 */
void MaterialFieldQuantity::DoComputeDebugCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int& cell_data_counter, MaterialName const material ) const {

   // extract the correct field type
   MaterialFieldType const field_type = quantity_data_.field_type_;

   // Check whether the given material is contained in the node
   if( node.ContainsMaterial( material ) ) {
      // local counter
      unsigned long long int local_counter = 0;
      //Loop through all components
      for( std::size_t component = 0; component < quantity_data_.field_indices_.size(); component++ ) {

         // Get the correct buffer and dimensionalization factor
         unsigned int const field_index                                 = quantity_data_.field_indices_[component];
         double const( &field_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( material ).GetFieldBuffer( field_type, field_index, buffer_type_ );
         double const dimensionalization_factor                         = unit_handler_.DimensionalizeValue( 1.0, MF::FieldUnit( field_type, field_index ) );

         // set the local counter on original cell_data_counter + the component
         local_counter = cell_data_counter + component;

         // Loop through all cells
         for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
            for( unsigned int j = 0; j < CC::TCY(); ++j ) {
               for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                  // add the dimensionalized value
                  cell_data[local_counter] = field_buffer[i][j][k] * dimensionalization_factor;

                  local_counter += quantity_data_.field_indices_.size();
               }
            }
         }
      }
   } else {
      //Loop through all components
      for( std::size_t component = 0; component < quantity_data_.field_indices_.size(); component++ ) {
         // Loop through all cells
         for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
            for( unsigned int j = 0; j < CC::TCY(); ++j ) {
               for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                  // add the dimensionalized value
                  cell_data[cell_data_counter++] = quantity_data_.debug_default_value_;
               }
            }
         }
      }
   }
}