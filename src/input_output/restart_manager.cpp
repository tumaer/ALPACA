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
#include "input_output/restart_manager.h"

#include <vector>

#include "block_definitions/block.h"
#include "user_specifications/compile_time_constants.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "enums/interface_tag_definition.h"
#include "communication/mpi_utilities.h"
#include "materials/material_definitions.h"
#include "utilities/buffer_operations.h"

namespace {
   /**
   * @brief Function that checks two parameter of the same type on consistency, i.e. that they are the same.
   * @param first_value,second_value The two values that should be checked on consistency.
   * @param error_name The name that is used in the error message, if both values are not the same.
   * @return Throws error message, if values are not the same.
   * @tparam T Type of the two parameter that are checked on consistency.
   */
   template<typename T>
   void CheckConsistency( T const& first_value, T const& second_value, std::string const& error_name ) {
      if( first_value != second_value ) {
         throw std::logic_error( error_name + " of input and restart file do not match (" + std::to_string( first_value ) + " vs " + std::to_string( second_value ) + ")!" );
      }
   }
}// namespace

/**
 * @brief Constructs the manager for writing and reading restart snapshots.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @param tree The tree to read and write the simulation data. It is only modified if the simulation is restored from a snapshot.
 * @param topology The topology to get and set information about the global structure of the simulation. It is only modified if the simulation is restored from a snapshot.
 * @param maximum_level The maximum level of the simulation.
 */
RestartManager::RestartManager( UnitHandler const& unit_handler, TopologyManager& topology_manager, Tree& tree, unsigned int const maximum_level ) :// Start initializer list
                                                                                                                                                     tree_( tree ),
                                                                                                                                                     topology_( topology_manager ),
                                                                                                                                                     logger_( LogWriter::Instance() ),
                                                                                                                                                     hdf5_manager_( Hdf5Manager::Instance() ),
                                                                                                                                                     maximum_level_( maximum_level ),
                                                                                                                                                     length_reference_( unit_handler.DimensionalizeValue( 1.0, UnitType::Length ) ),
                                                                                                                                                     density_reference_( unit_handler.DimensionalizeValue( 1.0, UnitType::Density ) ),
                                                                                                                                                     temperature_reference_( unit_handler.DimensionalizeValue( 1.0, UnitType::Temperature ) ),
                                                                                                                                                     velocity_reference_( unit_handler.DimensionalizeValue( 1.0, UnitType::Velocity ) ) {
   /** Empty besides initializer list */
}

/**
 * @brief Restores the simulation topology and tree from a restart file.
 * @param restore_filename The name of the file that should be used to restore the simulation.
 * @return The time at which the restart snapshot file was written.
 */
double RestartManager::RestoreSimulation( std::string const& restore_filename ) const {

   /** Open the hdf5 file for reading */
   hdf5_manager_.OpenFile( restore_filename, Hdf5Access::Read );

   /** Read metadata and check on consistency with current simulation setup */
   hdf5_manager_.OpenGroup( "simulation_data" );

   // Number of Dimensions
   unsigned int dimensions;
   hdf5_manager_.ReadAttributeScalar( "Dimensions", dimensions, H5T_NATIVE_UINT );
   CheckConsistency( dimensions, DTI( CC::DIM() ), "Compile Constants (Dimensions)" );

   // Number of Internal Cells
   unsigned int internal_cells;
   hdf5_manager_.ReadAttributeScalar( "InternalCells", internal_cells, H5T_NATIVE_UINT );
   CheckConsistency( internal_cells, CC::ICX(), "Compile Constants (Internal cells)" );

   // Number of Halo Cells
   unsigned int halo_size;
   hdf5_manager_.ReadAttributeScalar( "HaloSize", halo_size, H5T_NATIVE_UINT );
   CheckConsistency( halo_size, CC::HS(), "Compile Constants (Halo size)" );

   // Maximum level
   unsigned int maximum_level_file;
   hdf5_manager_.ReadAttributeScalar( "MaximumLevel", maximum_level_file, H5T_NATIVE_UINT );
   CheckConsistency( maximum_level_file, maximum_level_, "Maximum level" );

   // Length reference
   double length_reference;
   hdf5_manager_.ReadAttributeScalar( "LengthReference", length_reference, H5T_NATIVE_DOUBLE );
   CheckConsistency( length_reference, length_reference_, "Length reference" );

   // Velocity reference
   double velocity_reference;
   hdf5_manager_.ReadAttributeScalar( "VelocityReference", velocity_reference, H5T_NATIVE_DOUBLE );
   CheckConsistency( velocity_reference, velocity_reference_, "Velocity reference" );

   // Density reference
   double density_reference;
   hdf5_manager_.ReadAttributeScalar( "DensityReference", density_reference, H5T_NATIVE_DOUBLE );
   CheckConsistency( density_reference, density_reference_, "Density reference" );

   // Temperature reference
   double temperature_reference;
   hdf5_manager_.ReadAttributeScalar( "TemperatureReference", temperature_reference, H5T_NATIVE_DOUBLE );
   CheckConsistency( temperature_reference, temperature_reference_, "Temperature reference" );

   // Time (no consitency check required)
   double restart_time;
   hdf5_manager_.ReadAttributeScalar( "Time", restart_time, H5T_NATIVE_DOUBLE );

   /** Close all opened groups */
   hdf5_manager_.CloseGroup();

   /** Read the global information information of all nodes (required to restore the topology) */
   hdf5_manager_.OpenGroup( "node_info_data" );

   // Get the information for the number of nodes of the full topology and then read the nodes ids
   unsigned int const global_number_of_nodes = hdf5_manager_.GetDatasetExtent( "NodeIds" );
   std::vector<nid_t> node_ids( global_number_of_nodes );
   hdf5_manager_.ReadFullDataset( "NodeIds", node_ids.data(), H5T_NATIVE_ULLONG );

   // Get the information of all other required parameters (Consistency check on the given input data)
   std::vector<unsigned short> number_of_materials( global_number_of_nodes );
   hdf5_manager_.ReadFullDataset( "NumberOfMaterials", number_of_materials.data(), H5T_NATIVE_USHORT );

   std::vector<unsigned short> materials( std::accumulate( number_of_materials.begin(), number_of_materials.end(), 0 ) );
   hdf5_manager_.ReadFullDataset( "Materials", materials.data(), H5T_NATIVE_USHORT );

   std::vector<unsigned short> number_of_interface_blocks( global_number_of_nodes );
   hdf5_manager_.ReadFullDataset( "NumberOfInterfaceBlocks", number_of_interface_blocks.data(), H5T_NATIVE_USHORT );

   /** Close the open group */
   hdf5_manager_.CloseGroup();

   /** Restore the topology with the given data and obtain the local indices (of all global node ids) of the restored topology */
   auto const local_node_indices = topology_.RestoreTopology( node_ids, number_of_materials, materials );

   /** Now read all cell data from file and store it into the buffer */
   hdf5_manager_.OpenGroup( "node_cell_data" );

   // Define the size of the buffers that are read at once
   std::vector<hsize_t> const local_dimensions_conservatives( { 1, MF::ANOE(), CC::TCX(), CC::TCY(), CC::TCZ() } );
   std::vector<hsize_t> const local_dimensions_prime_states( { 1, MF::ANOP(), CC::TCX(), CC::TCY(), CC::TCZ() } );
   std::vector<hsize_t> const local_dimensions_single_buffer( { 1, CC::TCX(), CC::TCY(), CC::TCZ() } );

   // Open all datasets
   hdf5_manager_.OpenDatasetForReading( "Conservatives", local_dimensions_conservatives, H5T_NATIVE_DOUBLE );
   hdf5_manager_.OpenDatasetForReading( "PrimeStates", local_dimensions_prime_states, H5T_NATIVE_DOUBLE );
   hdf5_manager_.OpenDatasetForReading( "Levelset", local_dimensions_single_buffer, H5T_NATIVE_DOUBLE );
   hdf5_manager_.OpenDatasetForReading( "InterfaceTags", local_dimensions_single_buffer, H5T_NATIVE_CHAR );

   // Declare the buffers that are filled during reading plus other variables required during reading
   std::int8_t interface_tags[CC::TCX()][CC::TCY()][CC::TCZ()];
   double single_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
   std::vector<MaterialName> materials_of_node;

   // Loop through all local nodes and add the data
   for( auto const node_index : local_node_indices ) {

      // Compute the offset where the current material blocks start
      hsize_t const material_block_offset = std::accumulate( number_of_materials.begin(), number_of_materials.begin() + node_index, 0 );

      materials_of_node.clear();
      for( unsigned int material_index = 0; material_index < number_of_materials[node_index]; ++material_index ) {
         materials_of_node.push_back( ITM( materials[material_block_offset + material_index] ) );
      }

      // Declare the interface block as null_pointer
      std::unique_ptr<InterfaceBlock> interface_block = nullptr;
      if( number_of_interface_blocks[node_index] == 1 ) {
         hsize_t const reading_offset = std::accumulate( number_of_interface_blocks.begin(), number_of_interface_blocks.begin() + node_index, 0 );

         // read the levelset
         hdf5_manager_.ReadDataset( "Levelset", single_buffer, reading_offset );
         // Create an interface block with the read levelset values
         interface_block = std::make_unique<InterfaceBlock>( single_buffer );

         // read the interface tags
         hdf5_manager_.ReadDataset( "InterfaceTags", interface_tags, reading_offset );

      } else {
         /** Since it is already checked that the number of materials and number of interface blocks correspond to a correct arrangements,
          *  no further checks are required
          */
         // for a single-phase node the interface tags are uniform
         std::int8_t const uniform_tag = MaterialSignCapsule::SignOfMaterial( materials_of_node.front() ) * ITTI( IT::BulkPhase );
         BO::SetSingleBuffer( interface_tags, uniform_tag );
      }

      // Create the node with the material and interface data
      Node& new_node = tree_.CreateNode( node_ids[node_index], materials_of_node, interface_tags, std::move( interface_block ) );

      // Read the conservative and prime state data
      for( unsigned int material_index = 0; material_index < number_of_materials[node_index]; ++material_index ) {
         MaterialName const material = ITM( materials[material_block_offset + material_index] );
         Block& material_block       = new_node.GetPhaseByMaterial( material );

         hsize_t const reading_offset = material_block_offset + material_index;

         // Read the conservatives and primes
         hdf5_manager_.ReadDataset( "Conservatives", &material_block.GetRightHandSideBuffer(), reading_offset );
         hdf5_manager_.ReadDataset( "PrimeStates", &material_block.GetPrimeStateBuffer(), reading_offset );
      }
   }

   /** Close the file (closes also all groups and datasets that are open) */
   hdf5_manager_.CloseFile();

   // return the time of the restart file
   return restart_time;
}

/**
 * @brief Writes a restart file containing the topology and all relevant node data at the current time step.
 * @param timestep The current time step where the restart snapshot file is written.
 * @param filename_without_extension The filename of the restart file (without file extension).
 * @return The final name of the restart file that has been written.
 */
std::string RestartManager::WriteRestartFile( double const timestep, std::string const& filename_without_extension ) const {

   /** Prepare data that is required for the restart file (data that needs to be collected from different ranks, etc.) */
   unsigned int const my_rank = MpiUtilities::MyRankId();

   // node and material block data
   std::pair<unsigned int, unsigned int> const nodes_blocks_offset = topology_.NodeAndBlockOffsetOfRank( my_rank );
   std::pair<unsigned int, unsigned int> const nodes_blocks_global = topology_.NodeAndBlockCount();

   // Deduce variables for better readability
   unsigned int const global_number_of_nodes           = nodes_blocks_global.first;
   unsigned int const local_nodes_offset               = nodes_blocks_offset.first;
   unsigned int const global_number_of_material_blocks = nodes_blocks_global.second;
   unsigned int const local_material_blocks_offset     = nodes_blocks_offset.second;

   // interface block data (So far nodes can have only one levelset, so this way of counting is fine).
   unsigned int const local_number_of_interface_blocks = tree_.NodesWithLevelset().size();
   std::vector<unsigned int> number_of_interface_blocks_per_rank( MpiUtilities::NumberOfRanks() );
   MPI_Allgather( &local_number_of_interface_blocks, 1, MPI_UNSIGNED, number_of_interface_blocks_per_rank.data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD );
   unsigned int const local_interface_block_offset      = std::accumulate( number_of_interface_blocks_per_rank.begin(), number_of_interface_blocks_per_rank.begin() + my_rank, 0u );
   unsigned int const global_number_of_interface_blocks = std::accumulate( number_of_interface_blocks_per_rank.begin() + my_rank, number_of_interface_blocks_per_rank.end(), local_interface_block_offset );

   /** Define the dimensions for the different values and datasets that are written to the restart file */
   std::vector<hsize_t> const total_dimensions_conservatives( { global_number_of_material_blocks, MF::ANOE(), CC::TCX(), CC::TCY(), CC::TCZ() } );
   std::vector<hsize_t> const local_dimensions_conservatives( { 1, MF::ANOE(), CC::TCX(), CC::TCY(), CC::TCZ() } );

   std::vector<hsize_t> const total_dimensions_prime_states( { global_number_of_material_blocks, MF::ANOP(), CC::TCX(), CC::TCY(), CC::TCZ() } );
   std::vector<hsize_t> const local_dimensions_prime_states( { 1, MF::ANOP(), CC::TCX(), CC::TCY(), CC::TCZ() } );

   std::vector<hsize_t> const total_dimensions_single_buffer( { global_number_of_interface_blocks, CC::TCX(), CC::TCY(), CC::TCZ() } );
   std::vector<hsize_t> const local_dimensions_single_buffer( { 1, CC::TCX(), CC::TCY(), CC::TCZ() } );

   std::vector<hsize_t> const total_dimensions_node_scalar( { global_number_of_nodes } );
   std::vector<hsize_t> const local_dimensions_node_scalar( { 1 } );

   std::vector<hsize_t> const total_dimensions_block_scalar( { global_number_of_material_blocks } );
   std::vector<hsize_t> const local_dimensions_block_scalar( { 1 } );

   /** Open the hdf5 file where the data is written into */
   std::string const filename = filename_without_extension + ".h5";
   hdf5_manager_.OpenFile( filename );

   /** Open the group and datasets, where the basic information of the node are written into */
   hdf5_manager_.OpenGroup( "node_info_data" );
   hdf5_manager_.OpenDatasetForWriting( "NodeIds", total_dimensions_node_scalar, local_dimensions_node_scalar, local_nodes_offset, H5T_NATIVE_ULLONG );
   hdf5_manager_.OpenDatasetForWriting( "NumberOfMaterials", total_dimensions_node_scalar, local_dimensions_node_scalar, local_nodes_offset, H5T_NATIVE_USHORT );
   hdf5_manager_.OpenDatasetForWriting( "Materials", total_dimensions_block_scalar, local_dimensions_block_scalar, local_material_blocks_offset, H5T_NATIVE_USHORT );
   hdf5_manager_.OpenDatasetForWriting( "NumberOfInterfaceBlocks", total_dimensions_node_scalar, local_dimensions_node_scalar, local_nodes_offset, H5T_NATIVE_USHORT );

   /** Open the group and datasets, where the cell data information of the node are written into */
   hdf5_manager_.OpenGroup( "node_cell_data" );
   hdf5_manager_.OpenDatasetForWriting( "Conservatives", total_dimensions_conservatives, local_dimensions_conservatives, local_material_blocks_offset, H5T_NATIVE_DOUBLE );
   hdf5_manager_.OpenDatasetForWriting( "PrimeStates", total_dimensions_prime_states, local_dimensions_prime_states, local_material_blocks_offset, H5T_NATIVE_DOUBLE );
   hdf5_manager_.OpenDatasetForWriting( "Levelset", total_dimensions_single_buffer, local_dimensions_single_buffer, local_interface_block_offset, H5T_NATIVE_DOUBLE );
   hdf5_manager_.OpenDatasetForWriting( "InterfaceTags", total_dimensions_single_buffer, local_dimensions_single_buffer, local_interface_block_offset, H5T_NATIVE_CHAR );

   /** Write all node data to the file */
   for( auto const& level : tree_.FullNodeList() ) {
      for( auto const& [id, node] : level ) {
         /** Write general node info data */
         std::unordered_map<MaterialName, Block> const& phases( node.GetPhases() );
         unsigned short const number_of_materials        = phases.size();
         unsigned short const number_of_interface_blocks = node.HasLevelset() ? 1 : 0;

         // Write the info data
         hdf5_manager_.WriteDataset( "NodeIds", &id );
         hdf5_manager_.WriteDataset( "NumberOfInterfaceBlocks", &number_of_interface_blocks );
         hdf5_manager_.WriteDataset( "NumberOfMaterials", &number_of_materials );
         for( auto const& mat_block : phases ) {
            unsigned short const material_index = MTI( mat_block.first );
            hdf5_manager_.WriteDataset( "Materials", &material_index );
         }

         /** Write the actual cell data */
         // Write material/block data (conservatives and prime states)
         for( auto const& mat_block : phases ) {
            // Write conservatives and prime states
            hdf5_manager_.WriteDataset( "Conservatives", &mat_block.second.GetAverageBuffer() );
            hdf5_manager_.WriteDataset( "PrimeStates", &mat_block.second.GetPrimeStateBuffer() );
         }
         // write interface data (levelset and interface tags)
         if( node.HasLevelset() ) {
            hdf5_manager_.WriteDataset( "Levelset", node.GetInterfaceBlock().GetBaseBuffer( InterfaceDescription::Levelset ) );
            hdf5_manager_.WriteDataset( "InterfaceTags", node.GetInterfaceTags() );
         }
      }
   }

   /** Close the open groups (automatically closes all datasets) */
   hdf5_manager_.CloseGroup();

   /** Write metadata into the hdf5 file (all required to check consistency of restart file and new input file) */
   hdf5_manager_.OpenGroup( "simulation_data" );
   // Time
   hdf5_manager_.WriteAttributeScalar( "Time", timestep, H5T_NATIVE_DOUBLE );
   // Number of Dimensions
   hdf5_manager_.WriteAttributeScalar( "Dimensions", DTI( CC::DIM() ), H5T_NATIVE_UINT );
   // Number of Internal Cells
   hdf5_manager_.WriteAttributeScalar( "InternalCells", CC::ICX(), H5T_NATIVE_UINT );
   // Number of Halo Cells
   hdf5_manager_.WriteAttributeScalar( "HaloSize", CC::HS(), H5T_NATIVE_UINT );
   // Maximum level
   hdf5_manager_.WriteAttributeScalar( "MaximumLevel", maximum_level_, H5T_NATIVE_UINT );
   // Length reference
   hdf5_manager_.WriteAttributeScalar( "LengthReference", length_reference_, H5T_NATIVE_DOUBLE );
   // Velocity reference
   hdf5_manager_.WriteAttributeScalar( "VelocityReference", velocity_reference_, H5T_NATIVE_DOUBLE );
   // Density reference
   hdf5_manager_.WriteAttributeScalar( "DensityReference", density_reference_, H5T_NATIVE_DOUBLE );
   // Temperature reference
   hdf5_manager_.WriteAttributeScalar( "TemperatureReference", temperature_reference_, H5T_NATIVE_DOUBLE );

   /** Close the file (automatically closes all groups and datasets) */
   hdf5_manager_.CloseFile();

   /** Return the created filename and provide logging */
   logger_.LogMessage( "Restart file written at t = " + StringOperations::ToScientificNotationString( timestep, 9 ), true, true );
   return filename;
}
