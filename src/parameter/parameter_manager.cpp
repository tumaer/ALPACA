#include "parameter/parameter_manager.h"

#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief Standard constructor to create the ParameterManager object
 * @param material_manager Instance for handling models of material proeprties (such as viscosity)
 */
ParameterManager::ParameterManager( MaterialManager const& material_manager ) : material_manager_( material_manager ) {
   /** Empty besides initializer list */
}

/**
 * @brief Updates all parameters of a single node 
 * @param node The node under consideration (indirect return)
 */
void ParameterManager::UpdateParameter( Node & node ) const {

   // Obtain the cell size for the given node 
   double const cell_size = node.GetCellSize();

   // Update all parameters acting on a single material 
   for ( auto &phase : node.GetPhases() ) {
   
      // Obtain the sign of the material and the material 
      std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial( phase.first );
      Material const& material = material_manager_.GetMaterial( phase.first );

      // Start block to update the material property models 
      // 1. Shear viscosity
      if constexpr( CC::ViscosityIsActive() && CC::ShearViscosityModelActive() ) {
         // Call appropriate function depending on presence of an interface of this node 
         if ( node.HasLevelset() ) {
            material.GetShearViscosityModel().UpdateParameter( phase.second, cell_size, node.GetInterfaceTags(), material_sign );
         } else {
            material.GetShearViscosityModel().UpdateParameter( phase.second, cell_size );
         }

      }

      // 2. Thermal conductivity
      if constexpr( CC::HeatConductionActive() && CC::ThermalConductivityModelActive() ) {
         // Call appropriate function depending on presence of an interface of this node 
         if ( node.HasLevelset() ) {
            material.GetThermalConductivityModel().UpdateParameter( phase.second, cell_size, node.GetInterfaceTags(), material_sign );
         } else {
            material.GetThermalConductivityModel().UpdateParameter( phase.second, cell_size );
         }
      }
   }
}