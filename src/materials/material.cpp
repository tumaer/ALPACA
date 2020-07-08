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
#include "materials/material.h"

#include "materials/material_property_definitions.h"
#include "utilities/string_operations.h"

/**
 * @brief Sets up a material for the given input data.
 * @param equation_of_state Unique pointer to the equation of state that is used for the material (ownership transfer takes place).
 * @param dimensional_specific_heat_capacity Fixed dimensional value of the specific heat capacity.
 * @param dimensional_bulk_viscosity Fixed dimensional value of the bulk viscosity.
 * @param dimensional_shear_viscosity Fixed dimensional value of the shear viscosity.
 * @param dimensional_thermal_conductivity Fixed dimensional value of the thermal conductivity.
 * @param shear_viscosity_model Parameter model of the shear viscosity.
 * @param thermal_conductivity_model Parameter model of the thermal conductivity (ownership transfer takes place).
 * @param unit_handler Instance to provide (non-)dimensionalization of values (ownership transfer takes place).
 *
 * @note No default values are specified to ensure that everything is provided. Default values are set during the initialization.
 */
Material::Material( std::unique_ptr<EquationOfState const> equation_of_state,
                    double const dimensional_bulk_viscosity,
                    double const dimensional_shear_viscosity,
                    double const dimensional_thermal_conductivity,
                    double const dimensional_specific_heat_capacity,
                    std::unique_ptr<MaterialParameterModel const> shear_viscosity_model,
                    std::unique_ptr<MaterialParameterModel const> thermal_conductivity_model,
                    UnitHandler const& unit_handler ) :// Start initializer list
                                                        equation_of_state_( std::move( equation_of_state ) ),
                                                        bulk_viscosity_( unit_handler.NonDimensionalizeValue( dimensional_bulk_viscosity, UnitType::Viscosity ) ),
                                                        shear_viscosity_( unit_handler.NonDimensionalizeValue( dimensional_shear_viscosity, UnitType::Viscosity ) ),
                                                        thermal_conductivity_( unit_handler.NonDimensionalizeValue( dimensional_thermal_conductivity, UnitType::ThermalConductivity ) ),
                                                        specific_heat_capacity_( unit_handler.NonDimensionalizeValue( dimensional_specific_heat_capacity, { UnitType::Velocity, UnitType::Velocity }, { UnitType::Temperature } ) ),
                                                        shear_viscosity_model_( std::move( shear_viscosity_model ) ),
                                                        thermal_conductivity_model_( std::move( thermal_conductivity_model ) ) {
   /** Empty besides initializer list */
}

/**
 * @brief Move constructor.
 * @param material Material from which the data are moved.
 */
Material::Material( Material&& material ) :// Start initializer list
                                            equation_of_state_( std::move( material.equation_of_state_ ) ),
                                            bulk_viscosity_( material.bulk_viscosity_ ),
                                            shear_viscosity_( material.shear_viscosity_ ),
                                            thermal_conductivity_( material.thermal_conductivity_ ),
                                            specific_heat_capacity_( material.specific_heat_capacity_ ),
                                            shear_viscosity_model_( std::move( material.shear_viscosity_model_ ) ),
                                            thermal_conductivity_model_( std::move( material.thermal_conductivity_model_ ) ) {
   /** Empty besides initializer list */
}

/**
 * @brief Gives an instance to the equation of state for the given material.
 * @return equation of state.
 */
EquationOfState const& Material::GetEquationOfState() const {
   return *equation_of_state_;
}

/**
 * @brief Gives the shear and bulk viscosity (fixed values).
 * @return shear and bulk viscosity.
 */
std::vector<double> Material::GetShearAndBulkViscosity() const {
   return { shear_viscosity_, bulk_viscosity_ };
}

/**
 * @brief Gives the shear viscosity (fixed value).
 * @return shear viscosity.
 */
double Material::GetShearViscosity() const {
   return shear_viscosity_;
}

/**
 * @brief Gives the bulk viscosity (fixed value).
 * @return bulk viscosity.
 */
double Material::GetBulkViscosity() const {
   return bulk_viscosity_;
}

/**
 * @brief Gives the thermal conductivity (fixed value).
 * @return thermal conductivity.
 */
double Material::GetThermalConductivity() const {
   return thermal_conductivity_;
}

/**
 * @brief Gives the specific heat capacity (fixed value).
 * @return specific heat capacity.
 */
double Material::GetSpecificHeatCapacity() const {
   return specific_heat_capacity_;
}

/**
 * @brief Gives the model of the shear viscosity.
 * @return The instance of the shear viscosity model.
 */
MaterialParameterModel const& Material::GetShearViscosityModel() const {
   return *shear_viscosity_model_;
}

/**
 * @brief Gives the model of the thermal conductivity.
 * @return The instance of the thermal conductivity model.
 */
MaterialParameterModel const& Material::GetThermalConductivityModel() const {
   return *thermal_conductivity_model_;
}