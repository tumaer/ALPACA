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
#ifndef MATERIAL_FIELD_QUANTITIES_DEFINITIONS_H
#define MATERIAL_FIELD_QUANTITIES_DEFINITIONS_H

#include <iostream>
#include <vector>

#include "block_definitions/field_material_definitions.h"
#include "block_definitions/field_details.h"

/**
 * @brief The MaterialFieldOutputName enum gives all possible choices for which an output can be written (conservatives, primes states and paramteres)
 *
 * @note: In order to append a MaterialFieldQuantity to the ouput follow the steps:
 *        1. Add The field quantity to MaterialFieldQuantityName in the appropriate group (arbitrary name)
 *        2. Add for this name an entry with default values to the function DataOfMaterialFieldQuantity below
 *           ( for scalar quantities follow e.g. mass, for vectorial quantities follow e.g. velocity)
 *        3. Add an output constants expression in user_specifications/output_constants.h
 *        4. Add the constructor call in the "initialization_output_writer.cpp" function.
 */
enum class MaterialFieldQuantityName {
   // conservatives
   Mass,
   Momentum,
   Energy,
   // prime states
   Density,
   Temperature,
   Pressure,
   Velocity,
   // parameters
   ShearViscosity,
   ThermalConductivity
};

/**
 * @brief Returns the suffix for a given conservative buffer type.
 * @param buffer_type Conservative buffer type identifier.
 * @return suffix for the given identifier.
 */
inline std::string SuffixOfConservativeBufferType( ConservativeBufferType const buffer_type ) {
   switch( buffer_type ) {
      case ConservativeBufferType::Average: {
         return "";
      }
      case ConservativeBufferType::RightHandSide: {
         return "_rhs";
      }
      default: {//last possibility ConservativeBufferType::Initial :
         return "_initial";
      }
   }
}

/**
 * @brief Defines a struct storing all necessary information to provide correct information for each MaterialFieldQuantity.
 */
struct MaterialFieldQuantityData {
   // material field type for which the output is written
   MaterialFieldType field_type_;
   // indices of the field type which should be included in the output (more than one element generates vectorial output)
   std::vector<unsigned int> field_indices_;
   // default value for the debug output
   double debug_default_value_;

   /**
    * @brief Gives the dimensions of the quantity.
    * @return Dimensions of the quantity.
    */
   std::array<unsigned int, 2> GetDimensions() const {
      return { static_cast<unsigned int>( field_indices_.size() ), 1 };
   }

   /**
    * @brief Verifies the quantity data on validity to eliminate non-active fields.
    * @return True if any of the fields are active, otherwise false.
    * @note Manipulates the member data in case some fields are not active.
    */
   bool VerifyQuantity() {
      // Local field type for lambda capture
      MaterialFieldType const field_type = field_type_;
      // Store current size for warning handling
      size_t current_size = field_indices_.size();
      // Erase all elements that are not active
      field_indices_.erase( std::remove_if( field_indices_.begin(), field_indices_.end(),
                                            [&field_type]( double const& index ) { return !MF::IsFieldActive( field_type, index ); } ),
                            field_indices_.end() );
      // Throw warning if some of the fields have been rejected
      if( current_size != field_indices_.size() ) {
         std::cerr << "Output Warning! Some of the desired material fields are not active! Continue with reduced set!" << std::endl;
      }

      // If no field index is left anymore return false otherwise true
      return !field_indices_.empty();
   }
};

/***********************************************************************************************************************/
/*                              Add here new material field quantities                                                    */
/***********************************************************************************************************************/
/**
 * @brief Gives the appropriate data for a material field output quantity to generate the output.
 * @param output_quantity The material field output quantity identifier.
 * @return Quantity data struct for the given name.
 */
inline MaterialFieldQuantityData DataOfMaterialFieldQuantity( MaterialFieldQuantityName const output_quantity ) {
   // General declaration of the struct
   MaterialFieldQuantityData quantity_data;
   // Differ between all quantities
   switch( output_quantity ) {
      // conservatives
      case MaterialFieldQuantityName::Mass: {
         quantity_data.field_type_          = MaterialFieldType::Conservatives;
         quantity_data.field_indices_       = { ETI( Equation::Mass ) };
         quantity_data.debug_default_value_ = -1.0;
         return quantity_data;
      }
      case MaterialFieldQuantityName::Momentum: {
         quantity_data.field_type_    = MaterialFieldType::Conservatives;
         quantity_data.field_indices_ = { ETI( Equation::MomentumX )
#if DIMENSION != 1
                                                ,
                                          ETI( Equation::MomentumY )
#endif
#if DIMENSION == 3
                                                ,
                                          ETI( Equation::MomentumZ )
#endif
         };
         quantity_data.debug_default_value_ = -1.0;
         return quantity_data;
      }
      case MaterialFieldQuantityName::Energy: {
         quantity_data.field_type_          = MaterialFieldType::Conservatives;
         quantity_data.field_indices_       = { ETI( Equation::Energy ) };
         quantity_data.debug_default_value_ = -1.0;
         return quantity_data;
      }
      // prime states
      case MaterialFieldQuantityName::Density: {
         quantity_data.field_type_          = MaterialFieldType::PrimeStates;
         quantity_data.field_indices_       = { PTI( PrimeState::Density ) };
         quantity_data.debug_default_value_ = -1.0;
         return quantity_data;
      }
      case MaterialFieldQuantityName::Temperature: {
         quantity_data.field_type_          = MaterialFieldType::PrimeStates;
         quantity_data.field_indices_       = { PTI( PrimeState::Temperature ) };
         quantity_data.debug_default_value_ = -1.0;
         return quantity_data;
      }
      case MaterialFieldQuantityName::Pressure: {
         quantity_data.field_type_          = MaterialFieldType::PrimeStates;
         quantity_data.field_indices_       = { PTI( PrimeState::Pressure ) };
         quantity_data.debug_default_value_ = -1.0;
         return quantity_data;
      }
      case MaterialFieldQuantityName::Velocity: {
         quantity_data.field_type_    = MaterialFieldType::PrimeStates;
         quantity_data.field_indices_ = { PTI( PrimeState::VelocityX )
#if DIMENSION != 1
                                                ,
                                          PTI( PrimeState::VelocityY )
#endif
#if DIMENSION == 3
                                                ,
                                          PTI( PrimeState::VelocityZ )
#endif
         };
         quantity_data.debug_default_value_ = -1.0;
         return quantity_data;
      }
      // parameter
      case MaterialFieldQuantityName::ShearViscosity: {
         quantity_data.field_type_          = MaterialFieldType::Parameters;
         quantity_data.field_indices_       = { PTI( Parameter::ShearViscosity ) };
         quantity_data.debug_default_value_ = -1.0;
         return quantity_data;
      }
      case MaterialFieldQuantityName::ThermalConductivity: {
         quantity_data.field_type_          = MaterialFieldType::Parameters;
         quantity_data.field_indices_       = { PTI( Parameter::ThermalConductivity ) };
         quantity_data.debug_default_value_ = -1.0;
         return quantity_data;
      }
      // default if nothing matches
      default: {
         throw std::logic_error( "Material field output quantity is not known!" );
      }
   }
}

#endif// MATERIAL_FIELD_QUANTITIES_DEFINITIONS_H