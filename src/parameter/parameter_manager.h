#ifndef PARAMETER_MANAGER_H
#define PARAMETER_MANAGER_H

#include "parameter/material_parameter_model.h"
#include "materials/material_manager.h"
#include "topology/node.h"

/**
 * @brief The ParameterManager handles all updates regarding the parameter buffers that lie on a single material block. In general, it is not restricted to 
 *        models, such as viscosity or thermal conductivity models, that act on material properties. It can handle all models that act on a specific parameter
 *        buffer. 
 */
class ParameterManager {
   // Instance for handling parameter calculation of material data (e.g., viscosity)
   MaterialManager const& material_manager_;

public:
   ParameterManager() = delete;
   explicit ParameterManager( MaterialManager const& material_manager );
   ~ParameterManager() = default;
   ParameterManager( ParameterManager const& ) = delete;
   ParameterManager& operator=( ParameterManager const& ) = delete;
   ParameterManager( ParameterManager&& ) = delete;
   ParameterManager& operator=( ParameterManager&& ) = delete;

   void UpdateParameter( Node& node ) const;
};

#endif // PARAMETER_MANAGER_H