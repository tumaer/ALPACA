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
#ifndef MATERIAL_PROPERTY_DEFINITIONS_H
#define MATERIAL_PROPERTY_DEFINITIONS_H

#include <string>
#include <stdexcept>

#include "utilities/string_operations.h"

/**
 * @brief The MaterialPropertyName defines the identifiers to all different pure material and material pairing properties that exist.
 *        It is used to provide a defined mapping for the properties between input, transfer and output.
 * @note The identifier is also used for the mapping of specific property parameter models.
 */
enum class MaterialProperty {
   ShearViscosity,
   BulkViscosity,
   ThermalConductivity,
   SpecificHeatCapacity,
   SurfaceTensionCoefficient
};

/**
 * @brief Converts the material property to a proper string that can be used for input and output of the data.
 * @param property The material property identifier.
 * @return first_capitalized Flag to indicate if the first letter should be capitalized or not.
 *
 * @note Do not change the string obtained for the given property, since this is also used for mapping the property to a corresponding
 *       parameter model. In case change both.
 */
inline std::string MaterialPropertyToString( MaterialProperty const property, bool const first_capitalized = false ) {

   switch( property ) {
      case MaterialProperty::SpecificHeatCapacity : {
         return first_capitalized ? "SpecificHeatCapacity" : "specificHeatCapacity";
      }
      case MaterialProperty::ThermalConductivity : {
         return first_capitalized ? "ThermalConductivity" : "thermalConductivity";
      }
      case MaterialProperty::ShearViscosity : {
         return first_capitalized ? "ShearViscosity" : "shearViscosity";
      }
      case MaterialProperty::BulkViscosity : {
         return first_capitalized ? "BulkViscosity" : "bulkViscosity";
      }
      case MaterialProperty::SurfaceTensionCoefficient : {
         return first_capitalized ? "SurfaceTensionCoefficient" : "surfaceTensionCoefficient";
      }
      default : {
         // if nothing matches throw error
         throw std::logic_error( "Material property not known!" );
      }
   }
}

/**
 * @brief The MaterialPropertyModelName class defines all parameter models designed for computing material properties on other given material quantities
 *        such as prime states.
 *
 * @note The name MUST follow the syntax material property + model (exception is NotUsed).
 */
enum class MaterialPropertyModelName {
   // Default name if model is not used
   NotUsed,
   // viscosity models
      // constant model
      ShearViscosityConstant,
      // Shear rate models
      ShearViscosityPowerLaw,
      ShearViscosityCross,
      ShearViscosityCarreauYasuda,

      // Temperature models
      ShearViscositySutherland,

   // Thermal conductivity models
      // constant model
      ThermalConductivityConstant,

   // Surface tension coefficient models
      // constant model
      SurfaceTensionCoefficientConstant
};

/**
 * @brief Converts a string to its corresponding MaterialPropertyModelName.
 * @param property Material property for which the model is used (e.g. ShearViscosity, ThermalConductivity).
 * @param model_name Name of the model that is desired (e.g., Constant, Carreau).
 * @return Name of the material property model.
 */
inline MaterialPropertyModelName StringToMaterialPropertyModel( MaterialProperty const property, std::string const& model_name ) {
   // transform string to upper case without spaces
   std::string const name_upper_case( StringOperations::ToUpperCaseWithoutSpaces( MaterialPropertyToString( property ) + model_name ) );
   // Viscosity models
   // constant
   if( name_upper_case == "SHEARVISCOSITYCONSTANT" ) { return MaterialPropertyModelName::ShearViscosityConstant; }
   // shear rate models
   else if( name_upper_case == "SHEARVISCOSITYPOWERLAW" ) { return MaterialPropertyModelName::ShearViscosityPowerLaw; }
   else if( name_upper_case == "SHEARVISCOSITYCROSS" ) { return MaterialPropertyModelName::ShearViscosityCross; }
   else if( name_upper_case == "SHEARVISCOSITYCARREAUYASUDA" ) { return MaterialPropertyModelName::ShearViscosityCarreauYasuda; }
   // temperature models
   else if( name_upper_case == "SHEARVISCOSITYSUTHERLAND" ) { return MaterialPropertyModelName::ShearViscositySutherland; }

   // conductivity models
   // constant
   else if( name_upper_case == "THERMALCONDUCTIVITYCONSTANT" ) { return MaterialPropertyModelName::ThermalConductivityConstant; }

   // surface tension coefficient models
   // constant
   else if( name_upper_case == "SURFACETENSIONCOEFFICIENTCONSTANT" ) { return MaterialPropertyModelName::SurfaceTensionCoefficientConstant; }

   // default behavior if nothing is known
   else { throw std::logic_error( "Material property model is not known!" ); }
}

/**
 * @brief Converts a MaterialPropertyModelName into a corresponding string
 * @param name The identifier of the material property model
 * @return String to be used
 */
inline std::string MaterialPropertyModelToString( MaterialPropertyModelName const name ) {

   switch( name ) {
      // All constant models
      case MaterialPropertyModelName::ShearViscosityConstant :
      case MaterialPropertyModelName::ThermalConductivityConstant :
      case MaterialPropertyModelName::SurfaceTensionCoefficientConstant : {
         return "Constant";
      }
      // All specific models
      case MaterialPropertyModelName::ShearViscosityPowerLaw : {
         return "PowerLaw";
      }
      case MaterialPropertyModelName::ShearViscosityCross : {
         return "Cross";
      }
      case MaterialPropertyModelName::ShearViscosityCarreauYasuda : {
         return "CarreauYasuda";
      }
      case MaterialPropertyModelName::ShearViscositySutherland : {
         return "Sutherland";
      }
      // Default model not used
      default : {
         return "Not Used";
      }
   }
}

#endif // MATERIAL_PROPERTY_DEFINITIONS_H