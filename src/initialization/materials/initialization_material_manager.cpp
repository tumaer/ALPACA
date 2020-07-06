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
#include "initialization/materials/initialization_material_manager.h"

#include "initialization/materials/initialization_material.h"
#include "initialization/materials/initialization_material_pairing.h"

namespace Initialization {

   /**
    * @brief Initializes the complete set of materials with the given input reader.
    * @param material_reader Reader that provides access to the material data of the input file.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return The fully initialized set of materials present in the simulation.
    */
   std::vector<Material> InitializeMaterials( MaterialReader const& material_reader, UnitHandler const& unit_handler ) {
      // read the number of materials present in the simulation
      unsigned int const number_of_materials( material_reader.ReadNumberOfMaterials() );

      // Check if too much materials have been specified
      if( number_of_materials + 1 >= MTI( MaterialName::MaterialOutOfBounds ) ) {
         throw std::invalid_argument( "This number of materials is not currently not supported!" );
      }

      // In case of capillary forces at least tow materials must be defined
      if( CC::CapillaryForcesActive() && number_of_materials < 2 ) {
         throw std::invalid_argument( "For the simulation with Capillary Forces at least two materials must be specified!" );
      }

      // declare vector that is returned and reserve enough memory
      std::vector<Material> materials;
      materials.reserve( number_of_materials );
      // Initialize each material individually
      for( unsigned int mat_index = 0; mat_index < number_of_materials; mat_index++ ) {
         materials.push_back( InitializeMaterial( mat_index, material_reader, unit_handler ) );
      }

      // return the created vector
      return materials;
   }

   /**
    * @brief Initializes the complete set of material pairings with the given input reader.
    * @param material_reader Reader that provides access to the material data of the input file.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return The fully initialized set of material pairings present in the simulation.
    */
   std::vector<MaterialPairing> InitializeMaterialPairings( MaterialReader const& material_reader, UnitHandler const& unit_handler ) {
      // First read the number of all materials
      unsigned int const number_of_materials = material_reader.ReadNumberOfMaterials();
      // Check if too much materials have been specified
      if( number_of_materials + 1 >= MTI( MaterialName::MaterialOutOfBounds ) ) {
         throw std::invalid_argument( "This number of materials is currently not supported!" );
      }

      // Generate all combinations of indices
      std::vector<std::vector<unsigned int>> pairing_indices( GetMaterialPairingIndices( number_of_materials ) );

      // Check if too much material pairing combinations have been specified
      if( pairing_indices.size() + 1 >= MPTI( MaterialPairingName::MaterialPairingOutOfBounds ) ) {
         throw std::invalid_argument( "This number of material pairings is currently not supported!" );
      }

      // declare vector that is returned
      std::vector<MaterialPairing> material_pairings;
      // Definition and Initialization of the materialPairing vector depending on Single or Multiphase ( for later correct access )
      // NOTE: This if else can be removed if the general multimaterial approach is incorporated in the framework or the call to fill surface tension coefficient
      //       in capillary forces calculator is removed
      // Single material
      if( pairing_indices.size() == 0 ) {
         material_pairings.reserve( 1 );
         material_pairings.push_back( MaterialPairing() );
      }
      // Multimaterial
      else {
         // Reserve enough entries in the material pairing vector
         material_pairings.reserve( pairing_indices.size() );
         // Initialize each pairing individually
         for( auto const& indices : pairing_indices ) {
            material_pairings.push_back( InitializeMaterialPairing( indices, material_reader, unit_handler ) );
         }
      }

      // return the created vector
      return material_pairings;
   }


   /**
    * @brief Initializes the complete material manager class with the given input reader.
    * @param input_reader Reader that provides access to the full data of the input file.
    * @param unit_handler Instance to provide (non-)dimensionalization of values.
    * @return The fully initialized MaterialManager class.
    */
   MaterialManager InitializeMaterialManager( InputReader const& input_reader, UnitHandler const& unit_handler ) {

      // First create and then move for proper logging inside (correct order)
      std::vector<Material> materials( InitializeMaterials( input_reader.GetMaterialReader(), unit_handler ) );
      std::vector<MaterialPairing> material_pairings( InitializeMaterialPairings( input_reader.GetMaterialReader(), unit_handler ) );

      // Log a final empty line
      LogWriter & logger = LogWriter::Instance();
      logger.LogMessage( " " );

      return MaterialManager( std::move( materials ), std::move( material_pairings ) );
   }
}