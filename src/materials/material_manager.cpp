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
#include "materials/material_manager.h"

#include <stdexcept>
#include "utilities/string_operations.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
 * @brief Sets up a MaterialManager for the given pure materials and material pairings.
 * @param materials All already initialized materials.
 * @param material_pairing_data All already initialized material pairings.
 */
MaterialManager::MaterialManager( std::vector<Material> materials,
                                  std::vector<MaterialPairing> material_pairings ) :// Start initializer list
                                                                                     materials_( std::move( materials ) ),
                                                                                     material_pairings_( std::move( material_pairings ) ),
                                                                                     pairing_offset_( GenerateMaterialPairingOffset( materials_.size() ) ) {

   // Instantiation of material sign capsule for two-material computations
   // AB 20-03-30: This concept could be used everywhere if material manager is present in class and must be implemented in case of
   //              multi-material implementation
   MaterialSignCapsule( ITM( 0 ), ITM( materials_.size() - 1 ) );
}

/**
 * @brief Maps the pairing of two materials to the correct position of the vector using the material enum class index and a pairing offset.
 * @param first_material, second_material materials for which the pairing should be obtained.
 * @note This mapping ensures that the same result is ensured regardless of the order of both materials in the function call.
 */
unsigned int MaterialManager::MapPairingToIndex( MaterialName const first_material, MaterialName const second_material ) const {

   // NOTE: This part is required to have a fall back behavior in case capillary forces are activated and only one material is defined.
   if( materials_.size() == 1 ) {
      return 0;
   }

#ifndef PERFORMANCE
   if( first_material == second_material ) {
      throw std::runtime_error( "Error! Material pairing can only be obtained for two distinct materials!" );
   }
#endif

   // Obtain the indices of the two materials in the materials vector ( implicitly done due to the underlying index of MaterialName enum )
   auto material_A = MTI( first_material );
   auto material_B = MTI( second_material );

   // return the mapping index ( conversion to unsigned int is safe, since the max function always returns at least index 1 )
   return material_A + material_B + pairing_offset_[std::max( material_A, material_B ) - 1];
}

/**
 * @brief Returns the total number of materials contained in the simulation.
 * @return Number of materials.
 */
std::size_t MaterialManager::GetNumberOfMaterials() const {
   return materials_.size();
}

/**
 * @brief Returns all materials contained in the simulation.
 * @return Vector with all material names.
 */
std::vector<MaterialName> MaterialManager::GetMaterialNames() const {
   // Declare vector that is returned
   std::vector<MaterialName> material_names;
   // fill vector with all names
   material_names.reserve( materials_.size() );
   for( unsigned int material_index = 0; material_index < materials_.size(); material_index++ ) {
      material_names.push_back( ITM( material_index ) );
   }

   return material_names;
}

/**
 * @brief Gives an instance to the material with the given index.
 * @param index Index of the material that should be returned.
 * @return Instance of the material.
 * @note No sanity check is done here that the index really exists. Safety is ensured if the indices are
 *       obtained with the GetNumberOfMaterials() function.
 */
Material const& MaterialManager::GetMaterial( std::size_t const index ) const {

   return materials_[index];
}

/**
 * @brief Gives an instance of the material with the given identifier.
 * @param material The material identifier.
 * @return Instance of the material.
 * @note No sanity check is done here, since in general all present material names that are created should be
 *       checked by the MaterialName enum class itself.
 */
Material const& MaterialManager::GetMaterial( MaterialName const material ) const {

   return materials_[MTI( material )];
}

/**
 * @brief Gives an instance of a material pairing for two different materials.
 * @param first_material, second_material Unique identifiers of the materials of interest.
 * @return Instance of the material pairing.
 */
MaterialPairing const& MaterialManager::GetMaterialPairing( MaterialName const first_material, MaterialName const second_material ) const {
   return material_pairings_[MapPairingToIndex( first_material, second_material )];
}

/**
 * @brief Generates the offset values for each material to provide proper mapping of the material indices towards material pairing indices.
 * @param number_of_materials Number of materials to be considered.
 * @return Vector with the appropriate offset values for the given number of materials.
 */
std::vector<int> MaterialManager::GenerateMaterialPairingOffset( size_t const number_of_materials ) const {
   // declare vector
   std::vector<int> pairing_offset;
   // Definition and Initialization of the materialPairing vector depending on Single or Multiphase ( for later correct access )
   // NOTE: This if else can be removed if the general multimaterial approach is incorporated in the framework.
   // Single material
   if( number_of_materials == 1 ) {
      pairing_offset.reserve( 1 );
      pairing_offset.push_back( 0 );
   }
   // Multimaterial
   else {
      // Define the pairing offset for later mapping of the pairing to the appropriate positions in vector
      pairing_offset.reserve( number_of_materials - 1 );

      // The offset represents the increases of the index to access material_pairings_ vector resulting from an increase in material numbers
      // compared to the previous number of materials. The pairing_offset is then used together with the index of the materials enum class to access the
      // correct position in the material_pairings vector.
      // 2-materials: -1, 3-materials: -1, 4-materials: 0, 5-materials: 2; 6-fludis: 5, 7-materials: 9, ...
      int offset = -1;
      pairing_offset.push_back( offset );
      for( size_t n = 3; n <= number_of_materials; n++ ) {
         // increment offset
         offset += n - 3;
         // add offset at correct position for number of materials
         pairing_offset.push_back( offset );
      }
   }

   // return the vector
   return pairing_offset;
}

// MaterialSignCapsule does not have a cpp file for definition
MaterialName MaterialSignCapsule::positive_material_;
MaterialName MaterialSignCapsule::negative_material_;