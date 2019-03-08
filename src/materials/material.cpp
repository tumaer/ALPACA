#include "materials/material.h"

#include "materials/material_property_definitions.h"
#include "utilities/string_operations.h"

/**
 * @brief Sets up a material for the given input data 
 * @param equation_of_state Unique pointer to the equation of state that is used for the material (ownership transfer takes place)
 * @param specific_heat_capacity Fixed value of the specific heat capacity (dimensional)
 * @param bulk_viscosity Fixed value of the bulk viscosity (dimensional)
 * @param shear_viscosity Fixed value of the shear viscosity (dimensional)
 * @param thermal_conductivity Fixed value of the thermal conductivity (dimensional)
 * @param shear_viscosity_model Parameter model of the shear viscosity 
 * @param thermal_conductivity_model Parameter model of the thermal conductivity (ownership transfer takes place)
 * @param unit_handler Instance to provide (non-)dimensionalization of values (ownership transfer takes place)
 * 
 * @note No default values are specified to ensure that everything is provided. Default values are set during the initialization 
 */
Material::Material( std::unique_ptr<EquationOfState const> equation_of_state, 
                    double const bulk_viscosity,
                    double const shear_viscosity, 
                    double const thermal_conductivity, 
                    double const specific_heat_capacity, 
                    std::unique_ptr<MaterialParameterModel const> shear_viscosity_model,
                    std::unique_ptr<MaterialParameterModel const> thermal_conductivity_model,
                    UnitHandler const& unit_handler ) : 
   // Start initializer list
   equation_of_state_( std::move( equation_of_state ) ),
   bulk_viscosity_( unit_handler.NonDimensionalizeValue( bulk_viscosity, UnitType::Viscosity ) ),
   shear_viscosity_( unit_handler.NonDimensionalizeValue( shear_viscosity, UnitType::Viscosity ) ),
   thermal_conductivity_( unit_handler.NonDimensionalizeValue( thermal_conductivity, UnitType::ThermalConductivity ) ),
   specific_heat_capacity_( unit_handler.NonDimensionalizeValue( specific_heat_capacity, { UnitType::Velocity, UnitType::Velocity }, { UnitType::Temperature } ) ), 
   shear_viscosity_model_( std::move( shear_viscosity_model ) ),
   thermal_conductivity_model_( std::move( thermal_conductivity_model ) ) {
   /** Empty besides initializer list */
}

/**
 * @brief Move constructor
 * @param material Material from which the data are moved
 */
Material::Material( Material && material ) : 
   // Start initializer list
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
 * @brief Gives an instance to the equation of state for the given material
 * @return equation of state
 */
EquationOfState const& Material::GetEquationOfState() const {
   return *equation_of_state_;
}

/**
 * @brief Gives the shear and bulk viscosity (fixed values)
 * @return shear and bulk viscosity
 */
std::vector<double> Material::GetShearAndBulkViscosity() const {
   return { shear_viscosity_, bulk_viscosity_ };
}

/**
 * @brief Gives the shear viscosity (fixed value)
 * @return shear viscosity 
 */
double Material::GetShearViscosity() const {
   return shear_viscosity_;
}

/**
 * @brief Gives the bulk viscosity (fixed value)
 * @return bulk viscosity
 */
double Material::GetBulkViscosity() const {
   return bulk_viscosity_;
}

/**
 * @brief Gives the thermal conductivity (fixed value)
 * @return thermal conductivity
 */
double Material::GetThermalConductivity() const {
   return thermal_conductivity_;
}

/**
 * @brief Gives the specific heat capacity (fixed value)
 * @return specific heat capacity
 */
double Material::GetSpecificHeatCapacity() const {
   return specific_heat_capacity_;
}

/**
 * @brief Gives the model of the shear viscosity 
 * @return The instance of the shear viscosity model
 */
MaterialParameterModel const& Material::GetShearViscosityModel() const {
   return *shear_viscosity_model_;
}

/**
 * @brief Gives the model of the thermal conductivity 
 * @return The instance of the thermal conductivity model
 */
MaterialParameterModel const& Material::GetThermalConductivityModel() const {
   return *thermal_conductivity_model_;
}