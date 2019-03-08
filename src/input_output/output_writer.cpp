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
#include "input_output/output_writer.h"

#include "topology/id_information.h"
#include "communication/mpi_utilities.h"

#include "utilities/file_operations.h"
#include "utilities/string_operations.h"
#include "input_output/output_writer/xdmf_utilities.h"
#include "input_output/output_writer/output_definitions.h"

/**
 * @brief Creates an object to get the simulation data from the RAM to the hard disk.
 * @param standard_mesh_generator The already initialized standard mesh geenrator
 * @param debug_mesh_generator The already initialized debug mesh geenrator
 * @param interface_mesh_generator The already initialized interface mesh geenrator
 * @param material_output_quantities Vector holding all already initialized material output quantities that are written
 * @param interface_output_quantities Vector holding all already initialized interface output quantities that are written
 * @param number_of_materials Number of materials considered in the simulation
 *
 * @note for the pointer ownership transfer takes place
 */
OutputWriter::OutputWriter( std::unique_ptr<MeshGenerator const> standard_mesh_generator,
                            std::unique_ptr<MeshGenerator const> debug_mesh_generator,
                            std::unique_ptr<MeshGenerator const> interface_mesh_generator,
                            std::vector<std::unique_ptr<OutputQuantity const>> material_output_quantities,
                            std::vector<std::unique_ptr<OutputQuantity const>> interface_output_quantities,
                            unsigned int const number_of_materials ) :
   // Start initializer list
   number_of_materials_( number_of_materials ),
   hdf5_manager_( Hdf5Manager::Instance() ),
   standard_mesh_generator_( std::move( standard_mesh_generator ) ),
   debug_mesh_generator_( std::move( debug_mesh_generator ) ),
   interface_mesh_generator_( std::move( interface_mesh_generator ) ),
   material_output_quantities_( std::move( material_output_quantities ) ),
   interface_output_quantities_( std::move( interface_output_quantities ) ),
   material_quantities_dimension_map_( ComputeDimensionMap( material_output_quantities_ ) ),
   interface_quantities_dimension_map_( ComputeDimensionMap( interface_output_quantities_ ) ) {
   /** Empty besides initializer list */
}

/**
 * @brief Executes all necessary operations to initialize the xdmf time series file
 * @param time_series_filename_without_extension Name of the time series file (without extension)
 */
void OutputWriter::InitializeTimeSeriesFile( std::string const& time_series_filename_without_extension ) const {
   // Only do writing to file on 0 rank
   if( MpiUtilities::MyRankId() == 0 ) {
      // Get the header of the time series file
      std::string const xdmf_content = XdmfUtilities::HeaderInformation( "TimeSeries" );
      // Append content to the file of the time series
      FileOperations::WriteTextBasedFile( time_series_filename_without_extension + ".xdmf", xdmf_content );
   }
}

/**
 * @brief Executes all necessary operations to finalize the xdmf time series file
 * @param time_series_filename_without_extension Name of the time series file (without extension)
 */
void OutputWriter::FinalizeTimeSeriesFile( std::string const& time_series_filename_without_extension ) const {
   // Only do writing to file on 0 rank
   if( MpiUtilities::MyRankId() == 0 ) {
      // Get the footer of the time series file
      std::string const xdmf_content = XdmfUtilities::FooterInformation();
      // Append content to the file of the time series
      FileOperations::AppendToTextBasedFile( time_series_filename_without_extension + ".xdmf", xdmf_content );
   }
}

/**
 * @brief Triggers the output of the simulation results. Based on user Input the correct type of output is created.
 * @param output_type Output type identifier for which the output is written (standard, interface, debug)
 * @param output_time The time at which the simulation is currently at.
 * @param filename_without_extension Filename without extension where the output is written to
 * @param time_series_filename_without_extension Name of the time series file (without extension)
 */
void OutputWriter::WriteOutput( OutputType const output_type,
                                double const output_time,
                                std::string const& filename_without_extension,
                                std::string const& time_series_filename_without_extension ) const {

   // Define the full path filename of the hdf5 file
   std::string const hdf5_filename( filename_without_extension + ".h5" );
   // Obtain the correct mesh generator for the given output (ternary operator used to avoid new function declaration and allow constness)
   MeshGenerator const& mesh_generator = output_type == OutputType::Debug     ? *debug_mesh_generator_ :
                                         output_type == OutputType::Interface ? *interface_mesh_generator_ :
                                                                                *standard_mesh_generator_;
   // Call the writing function with the specific mesh_generator
   WriteHdf5File( output_time, hdf5_filename, mesh_generator, output_type );
   // Only write xdmf on rank 0 to avoid write conflicts
   if( MpiUtilities::MyRankId() == 0 ) {
      // Write the single time step xdmf file
      WriteXdmfTimeStepFile( output_time, hdf5_filename, mesh_generator, output_type );
      // Append the xdmf information to the time series file (ir required)
      if( !time_series_filename_without_extension.empty() ) {
         AppendToXdmfTimeSeriesFile( output_time, hdf5_filename, mesh_generator, time_series_filename_without_extension, output_type );
      }
   }
}

/**
 * @brief Writes a single time step xdmf file for a given time step
 * @param output_time Time at which the output is written
 * @param hdf5_filename Filename of the hdf5 file (absolute path), where the data information was written
 * @param mesh_generator The mesh generator instance which was used to write the output
 * @param output_type Type of the output to be written (standard, interface, debug)
 */
void OutputWriter::WriteXdmfTimeStepFile( double const output_time,
                                          std::string const& hdf5_filename,
                                          MeshGenerator const& mesh_generator,
                                          OutputType const output_type ) const {

   // Append header of the file
   std::string xdmf_content( XdmfUtilities::HeaderInformation( "TimeStep" ) );
   // Append core information of the xdmf file (topology, geometry and cell data)
   xdmf_content += XdmfSpatialDataInformation( output_time, FileOperations::RemoveFilePath( hdf5_filename ), mesh_generator, output_type );
   // Append the footer information
   xdmf_content += XdmfUtilities::FooterInformation();

   // Write the complete content to the xdmf file
   FileOperations::WriteTextBasedFile( FileOperations::ChangeFileExtension( hdf5_filename, ".xdmf" ), xdmf_content );
}

/**
 * @brief Append the data to the time series xdmf file
 * @param output_time Time at which the output is written
 * @param hdf5_filename Hdf5 file (absolute path), where the data information was written
 * @param mesh_generator The mesh generator instance which was used to write the output
 * @param output_type Type of the output to be written (standard, interface, debug)
 */
void OutputWriter::AppendToXdmfTimeSeriesFile( double const output_time,
                                               std::string const& hdf5_filename,
                                               MeshGenerator const& mesh_generator,
                                               std::string const& time_series_filename_without_extension,
                                               OutputType const output_type ) const {

   // Generate the core which should be appended to the file (use incrementing grid names)
   std::string const xdmf_content = XdmfSpatialDataInformation( output_time, FileOperations::RemoveFilePath( hdf5_filename ), mesh_generator, output_type );
   // Append content to the file of the time series
   FileOperations::AppendToTextBasedFile( time_series_filename_without_extension + ".xdmf", xdmf_content );
}

/**
 * @brief Writes the information contained in the xdmf file
 * @param output_time Time where the output is written
 * @param hdf5_short_filename short filename of the hdf5 file (without path information)
 * @param mesh_generator The mesh generator instance which was used to write the output
 * @param output_type Type of the output to be used (standard, interface, debug)
 * @return String of the Xdmf data that should be written
 */
std::string OutputWriter::XdmfSpatialDataInformation( double const output_time,
                                                      std::string const& hdf5_short_filename,
                                                      MeshGenerator const& mesh_generator,
                                                      OutputType const output_type ) const {

   // Declare the grid name with a prefix and given time used in the name of the hdf5 file
   std::string const spatial_data_name( "SpatialData_" + StringOperations::ToScientificNotationString( output_time, 6 ) );
   // Declare output string
   std::string xdmf_content;
   // Append the time information (empty line afterwards for visual separation)
   xdmf_content += XdmfUtilities::TimeDataItem( output_time );
   xdmf_content += '\n';
   // Append the Topology and domain information (empty line afterwards for visual separation)
   xdmf_content += mesh_generator.GetXdmfTopologyString( hdf5_short_filename, "mesh_topology" );
   xdmf_content += mesh_generator.GetXdmfGeometryString( hdf5_short_filename, "mesh_topology" );
   xdmf_content += '\n';
   // Append for all material quantities the appropriate information
   hsize_t const global_number_cells = mesh_generator.GetGlobalNumberOfCells();
   for( auto const& output_quantity : material_output_quantities_ ) {
      if( output_quantity->IsActive( output_type ) ) {
         // Differ between debug and other outputs
         if( output_type == OutputType::Debug ) {
            // Add material prefix for all quantities
            for( size_t material_index = 0; material_index < number_of_materials_; material_index++ ) {
               xdmf_content += output_quantity->GetXdmfAttributeString( hdf5_short_filename, "cell_data", global_number_cells, "material_" + std::to_string( material_index + 1 ) + "_");
            }
         } else {
            xdmf_content += output_quantity->GetXdmfAttributeString( hdf5_short_filename, "cell_data", global_number_cells );
         }
      }
   }
   // For interface quantities no distinction is required for the debug output
   for( auto const& output_quantity : interface_output_quantities_ ) {
      if( output_quantity->IsActive( output_type ) ) {
         xdmf_content += output_quantity->GetXdmfAttributeString( hdf5_short_filename, "cell_data", global_number_cells );
      }
   }
   // return the string including surrounding grid information
   return XdmfUtilities::SpatialDataInformation( spatial_data_name, xdmf_content );
}

/**
 * @brief Writes the data into the hdf5 file
 * @param hdf5_filename Filename of the hdf5 file (absolute path)
 * @param mesh_generator The mesh generator to be used for the output
 * @param output_type Type of the output that is considered (standard, interface, debug)
 */
void OutputWriter::WriteHdf5File( double const output_time, std::string const& hdf5_filename, MeshGenerator const& mesh_generator, OutputType const output_type ) const {

   /** Open the hdf5 file */
   hdf5_manager_.OpenFile( hdf5_filename );

   /** Write metadata into the hdf5 file (currently only time) */
   hdf5_manager_.OpenGroup( "metadata" );
   hdf5_manager_.WriteAttributeScalar( "time", output_time, H5T_NATIVE_DOUBLE );
   hdf5_manager_.CloseGroup();

   /** Write complete mesh topology information into the hdf5 file */
   // Open mesh topology group
   hdf5_manager_.OpenGroup( "mesh_topology" );

   /** Vertex IDs */
   // Define the hdf5 dataset properties for vertex ids
   std::vector<hsize_t> const global_dimensions_vertex_ids = mesh_generator.GetGlobalDimensionsOfVertexIDs();
   std::vector<hsize_t> const local_dimensions_vertex_ids = mesh_generator.GetLocalDimensionsOfVertexIDs();
   hsize_t const start_index_vertex_ids = mesh_generator.GetLocalVertexIDsStartIndex();

   // Declare the vector where the vertex IDs are written into. The vector is not initialized with a certain size since done inside of the mesh generator
   std::vector<unsigned long long int> vertex_ids;
   // Compute the vertex coordinates
   mesh_generator.ComputeVertexIDs( vertex_ids );
   // Reserve storage for the vertex IDs dataset and define hyperslab positions
   hdf5_manager_.ReserveDataspace( "VertexIDs", global_dimensions_vertex_ids, local_dimensions_vertex_ids, start_index_vertex_ids, H5T_NATIVE_ULLONG );
   // Write data to the dataset
   hdf5_manager_.WriteDatasetToDataspace( "VertexIDs", mesh_generator.GetVertexIDsName(), vertex_ids.data() );
   // Release the memory for the dataset
   hdf5_manager_.CloseDataset( "VertexIDs" );

   /** Vertex coordinates */
   // Define the hdf5 dataset properties for vertex coordinates
   std::vector<hsize_t> const global_dimensions_vertex_coordinates = mesh_generator.GetGlobalDimensionsOfVertexCoordinates();
   std::vector<hsize_t> const local_dimensions_vertex_coordinates = mesh_generator.GetLocalDimensionsOfVertexCoordinates();
   hsize_t const start_index_vertex_coordinates = mesh_generator.GetLocalVertexCoordinatesStartIndex();

   // Declare the vector where the vertex coordinates are written into. The vector is not initialized with a certain size since done inside of the mesh generator
   std::vector<double> vertex_coordinates;
   // Compute the vertex coordinates
   mesh_generator.ComputeVertexCoordinates( vertex_coordinates );
   // Reserve storage for the vertex coordinates dataset and define hyperslab positions
   hdf5_manager_.ReserveDataspace( "VertexCoordinates", global_dimensions_vertex_coordinates, local_dimensions_vertex_coordinates, start_index_vertex_coordinates, H5T_NATIVE_DOUBLE );
   // Write data to the dataset
   hdf5_manager_.WriteDatasetToDataspace( "VertexCoordinates", mesh_generator.GetVertexCoordinatesName(), vertex_coordinates.data() );

   /** 2.4 Close domain group (automatically closes open datasets) */
   hdf5_manager_.CloseGroup();

   /** Write cell fields into the hdf5 file */
   /** Open group where cell data is written to*/
   hdf5_manager_.OpenGroup( "cell_data" );

   /** Define parameters used for all cell fields */
   // Local nodes that are written to the hdf5 file by the current rank
   std::vector<std::reference_wrapper<Node const>> const local_nodes = mesh_generator.GetLocalNodes();

   // Define the hdf5 dataset properties for cell data (same for all quantities)
   hsize_t const global_number_of_cells = mesh_generator.GetGlobalNumberOfCells();
   hsize_t const local_number_of_cells = mesh_generator.GetLocalNumberOfCells();
   hsize_t const local_cells_start_index = mesh_generator.GetLocalCellsStartIndex();

   // Data vector for each quantity (not initialized here. Size is assigned in the dimension loop)
   std::vector<double> cell_data;

   // Loop through all different material quantities dimensions
   for( auto const& [dimensions, quantity_indices] : material_quantities_dimension_map_ ) {

      // Data vector for each cell data (done here to avoid multiple calls to allocate and deallocate memory in each quantity)
      cell_data.resize( local_number_of_cells * dimensions[0] * dimensions[1] );

      // Reserve the dataset dataspace for all quantities of this dimension
      std::vector<hsize_t> const global_dimensions( { global_number_of_cells, dimensions[0], dimensions[1] } );
      std::vector<hsize_t> const local_dimensions( { local_number_of_cells, dimensions[0], dimensions[1] } );
      hdf5_manager_.ReserveDataspace( "BlockCellData", global_dimensions, local_dimensions, local_cells_start_index, H5T_NATIVE_DOUBLE );

      // Loop through all quantities with the given dimension
      for( auto const& quantity_index : quantity_indices ) {
         // Obtain the correct output quantity
         auto const& output_quantity = material_output_quantities_[quantity_index];
         // Check if the quantity is active for the given output type
         if( output_quantity->IsActive( output_type ) ) {
            // Change behavior of file writing dependent on output type
            if( output_type == OutputType::Debug ) {
               // Loop through all materials given in the current simulation
               for( size_t material_index = 0; material_index < number_of_materials_; material_index++ ) {
                  // Append the quantity to the hdf5 file
                  output_quantity->ComputeDebugCellData( local_nodes, cell_data, ITM( material_index ) );
                  // Open the dataset with the appropriate name
                  hdf5_manager_.WriteDatasetToDataspace( "BlockCellData", "material_" + std::to_string( material_index + 1 ) + "_" + output_quantity->GetName(), cell_data.data() );
               }
            } else {
               // Compute the cell data of the given quantity
               output_quantity->ComputeCellData( local_nodes, cell_data );
               // Write data to hdf5 file
               hdf5_manager_.WriteDatasetToDataspace( "BlockCellData", output_quantity->GetName(), cell_data.data() );
            }
         }
      }

      // Release the dataset dataspace
      hdf5_manager_.CloseDataset( "BlockCellData" );
   }

   // Loop through all different material quantities dimensions
   for( auto const& [dimensions, quantity_indices] : interface_quantities_dimension_map_ ) {

      // Data vector for each cell data (done here to avoid multiple calls to allocate and deallocate memory in each quantity)
      cell_data.resize( local_number_of_cells * dimensions[0] * dimensions[1] );

      // Reserve the dataset dataspace for all quantities of this dimension
      std::vector<hsize_t> const global_dimensions( { global_number_of_cells, dimensions[0], dimensions[1] } );
      std::vector<hsize_t> const local_dimensions( { local_number_of_cells, dimensions[0], dimensions[1] } );
      hdf5_manager_.ReserveDataspace( "InterfaceBlockCellData", global_dimensions, local_dimensions, local_cells_start_index, H5T_NATIVE_DOUBLE );

      // Loop through all quantities with the given dimension
      for( auto const& quantity_index : quantity_indices ) {
         // Obtain the correct output quantity
         auto const& output_quantity = interface_output_quantities_[quantity_index];
         // Check if the quantity is active for the given output type
         if( output_quantity->IsActive( output_type ) ) {
            // Compute the cell data depending on the given output type
            if( output_type == OutputType::Debug ) {
               output_quantity->ComputeDebugCellData( local_nodes, cell_data );
            } else {
               output_quantity->ComputeCellData( local_nodes, cell_data );
            }

            // Open the dataset with the appropriate name
            hdf5_manager_.WriteDatasetToDataspace( "InterfaceBlockCellData", output_quantity->GetName(), cell_data.data() );
         }
      }

      // Release the dataset dataspace
      hdf5_manager_.CloseDataset( "InterfaceBlockCellData" );
   }

   /** Close group */
   hdf5_manager_.CloseGroup();

   /** Closing the last HDF Ressources and write xdmf file */
   hdf5_manager_.CloseFile();
}

/**
 * @brief Computes a map that assigns for all different quantity dimensions the appropriate quantity indices in the output_quantity vector
 * @param quantities Vector with all quantities for which the map should be computed
 * @return Map with key: Array with both quantity dimensions (row, column); value: Vector with all indices in the quantities vector for the given dimension
 */
std::map<std::array<unsigned int,2>,std::vector<unsigned int>> OutputWriter::ComputeDimensionMap( std::vector<std::unique_ptr<OutputQuantity const>> const& quantities ) const {
   // Map that is returned
   std::map<std::array<unsigned int,2>,std::vector<unsigned int>> dimension_map;
   // quantity counter to get the correct index for each quantity
   unsigned int quantity_counter = 0;
   // Loop through all quantities
   for( auto const& quantity : quantities ) {
      // Create a new entry in the map for the given dimension or add the given index to existing
      dimension_map[ quantity->GetDimensions() ].push_back( quantity_counter );
      quantity_counter++;
   }

   return dimension_map;
}