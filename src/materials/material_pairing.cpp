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
#include "materials/material_pairing.h"

#include "materials/material_property_definitions.h"
#include "utilities/string_operations.h"

/**
 * @brief Sets up a material pairing for the given input data.
 * @param dimensional_surface_tension_coefficient Fixed dimensional value of the surface tension coefficient (ownership takes place).
 * @param surface_tension_coefficient_model The model used for the surface tension coefficient.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 *
 * @note No default values are specified to ensure that everything is provided. Default values are set during the initialization.
 */
MaterialPairing::MaterialPairing( double const dimensional_surface_tension_coefficient,
                                  std::unique_ptr<InterfaceParameterModel const> surface_tension_coefficient_model,
                                  UnitHandler const& unit_handler ) :// Start initializer list
                                                                      surface_tension_coefficient_( unit_handler.NonDimensionalizeValue( dimensional_surface_tension_coefficient, UnitType::SurfaceTensionCoefficient ) ),
                                                                      surface_tension_coefficient_model_( std::move( surface_tension_coefficient_model ) ) {
   /** Empty besides initializer list */
}

/**
 * @brief Sets up a material pairing with no input data
 * @note Function can be removed when the surface tension coefficient hard coding is removed from the capillary pressure calculator.
 */
MaterialPairing::MaterialPairing() : surface_tension_coefficient_( 0.0 ),
                                     surface_tension_coefficient_model_( nullptr ) {
   /** Empty besides initializer list */
}

/**
 * @brief Move constructor.
 * @param pairing Material pairing from which the data are moved.
 */
MaterialPairing::MaterialPairing( MaterialPairing&& pairing ) :// Start initializer list
                                                                surface_tension_coefficient_( pairing.surface_tension_coefficient_ ),
                                                                surface_tension_coefficient_model_( std::move( pairing.surface_tension_coefficient_model_ ) ) {
   /** Empty besides initializer list */
}

/**
 * @brief Gives the surface tension coefficient  (fixed value)
 * @return surface tension coefficient
 */
double MaterialPairing::GetSurfaceTensionCoefficient() const {
   return surface_tension_coefficient_;
}

/**
 * @brief Gives the model of the surface tension coefficient.
 * @return The instance of the surface tension coefficient model.
 */
InterfaceParameterModel const& MaterialPairing::GetSurfaceTensionCoefficientModel() const {
   return *surface_tension_coefficient_model_;
}
