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
#include "input_output/hdf5/hdf5_manager.h"

#include <stdexcept>

/**
 * @brief Default constructor (private).
 * @note can only be used to create a singleton hdf5 writer.
 */
Hdf5Manager::Hdf5Manager() {
   /** Empty constructor */
}

/**
 * @brief Destructor checks whether a file is still open and closes it.
 */
Hdf5Manager::~Hdf5Manager() {
   if( file_.is_open_ ) {
      CloseFile();
   }
}

/**
 * @brief Opens a file to write/read hdf5 content into/from.
 * @param filename Name of the file.
 * @param file_access The mode used to open the hdf5 file.
 */
void Hdf5Manager::OpenFile( std::string const& filename, Hdf5Access const file_access ) {
#ifndef PERFORMANCE
   // Check whether a file is already opened
   if( file_.is_open_ ) {
      throw std::runtime_error( "Error Opening HDF5 file. Cannot open two files at the same time!" );
   }
#endif
   // Set the flag for subsequent calls
   file_.access_type_ = file_access;

   // instantiates the file_properties
   file_.properties_ = H5Pcreate( H5P_FILE_ACCESS );
   H5Pset_fapl_mpio( file_.properties_, MPI_COMM_WORLD, MPI_INFO_NULL );
   // Opens the file
   file_.id_ = file_.access_type_ == Hdf5Access::Read ? H5Fopen( filename.c_str(), H5F_ACC_RDONLY, file_.properties_ ) : H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_.properties_ );
   // Set flag that file is open
   file_.is_open_ = true;
}

/**
 * @brief Opens a group to access hdf5 content from that group.
 * @param group_name Name of the group where the data is accessed from/to.
 *
 * @note After opening a group this group is always set to the active group. In case, different groups need to be accessed after opening, the
 *       function ActivateGroup can be used.
 */
void Hdf5Manager::OpenGroup( std::string const& group_name ) {
#ifndef PERFORMANCE
   // Check whether a file is already opened where the group could be appended
   if( !file_.is_open_ ) {
      throw std::runtime_error( "Error opening group. No file is open!" );
   }
   // Check if the group already exists
   if( groups_.find( group_name ) != groups_.end() ) {
      throw std::runtime_error( "Error opening group. The group already exists!" );
   }
#endif
   // Create locally the struct holding all group information
   Hdf5Group group;
   // open the group and properties used for accessing/writing datasets of/into the group
   group.id_         = file_.access_type_ == Hdf5Access::Read ? H5Gopen2( file_.id_, group_name.c_str(), H5P_DEFAULT ) : H5Gcreate2( file_.id_, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
   group.properties_ = H5Pcreate( H5P_DATASET_XFER );
   H5Pset_dxpl_mpio( group.properties_, H5FD_MPIO_INDEPENDENT );
   // Add the new group to the map
   groups_[group_name] = group;
   // Set the active group to the current opened group
   active_group_ = group_name;
}

/**
 * @brief Activates the group with given tag to open datasets/dataspaces into it.
 * @param name name of the group that should be activated.
 */
void Hdf5Manager::ActivateGroup( std::string const& name ) {
#ifndef PERFORMANCE
   // Check if the group exists
   if( groups_.find( name ) != groups_.end() ) {
      throw std::runtime_error( "Before activating a group, it must be opened!" );
   }
#endif
   active_group_ = name;
}

/**
 * @brief Gives the extent of an already existing dataset. The extent is defined by all number of points the dataset holds if it would be written into a single
 *        contiguous buffer.
 * @return The number of elements in the dataset.
 */
hsize_t Hdf5Manager::GetDatasetExtent( std::string const& dataset_name ) const {
#ifndef PERFORMANCE
   // Check if a group is active
   if( active_group_.empty() ) {
      throw std::runtime_error( "Before getting the extent of a dataset, a group must be opened and activated!" );
   }
#endif
   // Get the correct group info
   Hdf5Group const& group = groups_.at( active_group_ );

   /** Open the dataset for reading information from it */
   hid_t const dataset_id      = H5Dopen2( group.id_, dataset_name.c_str(), H5P_DEFAULT );
   hid_t const local_hyperslab = H5Dget_space( dataset_id );
   // read the data from file and store it in the contiguous buffer (cannot be returned directly, since data needs to be closed afterwards)
   hsize_t const dataset_extent = H5Sget_simple_extent_npoints( local_hyperslab );
   // close the hyperslab and dataset
   H5Sclose( local_hyperslab );
   H5Dclose( dataset_id );

   return dataset_extent;
}

/**
 * @brief Opens a dataset for writing data into it. It is part of a two-phase writing procedure:
 *           1. OpenDatasetForWriting (allocates the memory for the dataset to be written)
 *           2. WriteDataset (selects the hyperslab position and writes data into the dataset)
 *        Therefore, this function MUST be called before writing to a dataset.
 * @param dataset_name Name of the dataset that should be linked into the group/file.
 * @param dataspace_total_dimensions One-dimensional vector specifying the global dimensions of the dataspace (number of values + dimensions of each value)
 *        (e.g., 1000 cells where each is a 3D vector => {1000, 3, 1}, 1000 cells where each is a matrix nxm => {1000,n,m}).
 * @param dataset_local_dimensions The actual dimensions that are written to the file in one chunk
 *        (e.g. Global dimensions: {1000,3,1} where each vector is written once => local_dimensions = {1,3,1}).
 * @param local_elements_start_index The start index, where the current rank starts to write the cell values.
 * @param datatype_id Hdf5 identifier which datatype should be used for the dataset (e.g., H5T_NATIVE_ULLONG, H5_NATIVE_DOUBLE).
 *
 * @note The dataset is always appended to the current active group. Use ActivateGroup function to control.
 * @note The two-way writing procedure should be used, where already present (contiguous) data can be piped into the file without any further modifications. For each dataset, this function
 *       must be called. Two different datasets cannot be written in the same dataset.
 */
void Hdf5Manager::OpenDatasetForWriting( std::string const& dataset_name, std::vector<hsize_t> const& dataspace_total_dimensions, std::vector<hsize_t> const& dataset_local_dimensions, hsize_t const local_elements_start_index, hid_t const datatype_id ) {
#ifndef PERFORMANCE
   // Check if the active group flag is set
   if( active_group_.empty() ) {
      throw std::runtime_error( "Before opening a dataset for writing, a group must be active!" );
   }
   if( datasets_.find( dataset_name ) != datasets_.end() ) {
      throw std::runtime_error( "Error opening dataset for writing. The dataset already exists!" );
   }
#endif
   // Declare the dataset info struct to store all required information
   Hdf5Dataset dataset;

   /** Link the current active group to the dataset */
   dataset.group_name_    = active_group_;
   Hdf5Group const& group = groups_[dataset.group_name_];

   /** Copy and define the dimensional information about the dataset (done to have all information collectively stored in the struct) */
   dataset.total_dimensions_ = dataspace_total_dimensions;
   dataset.local_dimensions_ = dataset_local_dimensions;
   // Only the first start index is defined (represents the node, block, buffer position). All other are set to zero, since always a full local buffer is written.
   dataset.start_indices_         = std::vector<hsize_t>( dataset.local_dimensions_.size(), 0 );
   dataset.start_indices_.front() = local_elements_start_index;
   dataset.datatype_              = datatype_id;

   /** Define the chunk that is written at once and the proeprties used to create the dataset */
   // Define the chunk that is created at once (The same as local dimensions, but only one node/block/buffer is written once )
   dataset.chunk_             = dataset.total_dimensions_;
   dataset.chunk_.front()     = 1;
   dataset.properties_create_ = H5Pcreate( H5P_DATASET_CREATE );
   H5Pset_chunk( dataset.properties_create_, dataset.chunk_.size(), dataset.chunk_.data() );

   // Define the number of elements of the dataset that are written at once (per dimension). Always 1 per dimension.
   dataset.count_ = std::vector( dataset.total_dimensions_.size(), 1ull );

   /** Open the dataset for writing and select the entire appropriate hyperslab */
   // Create the dataspace for the full (global, all ranks) data (NULL marks that the dataspace is fixed to the given dimensions and cannot grow during writing)
   dataset.dataspace_id_ = H5Screate_simple( dataset.total_dimensions_.size(), dataset.total_dimensions_.data(), NULL );
   // Creates a new dataset and links it into the file/group
   dataset.dataset_id_ = H5Dcreate2( group.id_, dataset_name.c_str(), dataset.datatype_, dataset.dataspace_id_, H5P_DEFAULT, dataset.properties_create_, H5P_DEFAULT );
   // Create the local memory and hyperslab space (NULL marks that dataset cannot grow until it is closed)
   dataset.local_memory_space_ = H5Screate_simple( dataset.local_dimensions_.size(), dataset.local_dimensions_.data(), NULL );
   dataset.local_hyperslab_    = H5Dget_space( dataset.dataset_id_ );

   /** Add the new dataset to the map */
   datasets_[dataset_name] = dataset;
}

/**
 * @brief Opens a dataset for reading data from it. It is part of a two-phase reading procedure:
 *           1. OpenDatasetForReading (allocates the memory for the dataset to be read)
 *           2. ReadDataset (selects the hyperslab position and read data from the dataset)
 *        Therefore, this function MUST be called before reading from a dataset.
 * @param dataset_name Name of the dataset that should be read.
 * @param dataset_local_dimensions The actual dimensions that are read from the file in one chunk
 *        (e.g. Global dimensions: {1000,3,1} where each vector is read once => local_dimensions = {1,3,1}).
 * @param datatype_id Hdf5 identifier which datatype should be used for the dataset (e.g., H5T_NATIVE_ULLONG, H5_NATIVE_DOUBLE).
 *
 * @note The dataset is always appended to the current active group. Use ActivateGroup function to control.
 */
void Hdf5Manager::OpenDatasetForReading( std::string const& dataset_name, std::vector<hsize_t> const& dataset_local_dimensions, hid_t const datatype_id ) {
#ifndef PERFORMANCE
   // Check if the active group flag is set
   if( active_group_.empty() ) {
      throw std::runtime_error( "Before opening a dataset for reading, a group must be active!" );
   }
   if( datasets_.find( dataset_name ) != datasets_.end() ) {
      throw std::runtime_error( "Error opening dataset for reading. The dataset already exists!" );
   }
#endif
   // Declare the dataset info struct to store all required information
   Hdf5Dataset dataset;

   /** Link the current active group to the dataset */
   dataset.group_name_    = active_group_;
   Hdf5Group const& group = groups_[dataset.group_name_];

   /** Copy and define the dimensional information about the dataset */
   // Defines the local dimensions that are written to the global dataset dataspace
   dataset.local_dimensions_ = dataset_local_dimensions;
   // Defines the start indices for the local dimensions in the global dimensions matrix
   dataset.start_indices_ = std::vector<hsize_t>( dataset.local_dimensions_.size(), 0 );
   // Save the dataset data type
   dataset.datatype_ = datatype_id;

   /** Create the properties that are used for the reading procedure */
   // Define the count for the dataset
   dataset.count_ = std::vector( dataset.local_dimensions_.size(), 1ull );

   /** Open the dataset for reading and select the entire local hyperslab */
   // Creates a new dataset and links it into the file/group
   dataset.dataset_id_ = H5Dopen2( group.id_, dataset_name.c_str(), H5P_DEFAULT );
   // Create the local memory and hyperslab space required for writing the data (NULL marks that dataset cannot grow until it is closed)
   dataset.local_memory_space_ = H5Screate_simple( dataset.local_dimensions_.size(), dataset.local_dimensions_.data(), NULL );
   dataset.local_hyperslab_    = H5Dget_space( dataset.dataset_id_ );

   /** Add the new dataset to the map */
   datasets_[dataset_name] = dataset;
}

/**
 * @brief Reserves the dataspace for a dataset for the writing process. It is part of a two-phase writing procedure:
 *           1. ReserveDataspace (allocates the memory for the dataspace to be written)
 *           2. WriteDatasetToDataspace (creates the dataset, selects the hyperslab position and writes data to the dataset)
 *        Therefore, this function MUST be called before writing to a dataset.
 * @param dataspace_name Name of the dataspace that should be used.
 * @param dataspace_total_dimensions One-dimensional vector specifying the global dimensions of the dataspace (number of values + dimensions of each value)
 *        (e.g., 1000 cells where each is a 3D vector => {1000, 3, 1}, 1000 cells where each is a matrix nxm => {1000,n,m}).
 * @param dataset_local_dimensions The actual dimensions that are written to the file in one chunk.
 *        (e.g. Global dimensions: {1000,3,1} where each vector is written once => local_dimensions = {1,3,1}).
 * @param local_elements_start_index The start index, where the current rank starts to write the cell values.
 * @param datatype_id Hdf5 identifier which datatype should be used for the dataset (e.g., H5T_NATIVE_ULLONG, H5_NATIVE_DOUBLE).
 *
 * @note The dataset is always appended to the current active group. Use ActivateGroup function to control.
 * @note The two-way writing procedure should be used where different datasets have the same extent and only the name of the datasets and the values differ. Hence, different datasets
 *       can be written with the same dataspace.
 */
void Hdf5Manager::ReserveDataspace( std::string const& dataspace_name, std::vector<hsize_t> const& dataspace_total_dimensions, std::vector<hsize_t> const& dataset_local_dimensions, hsize_t const local_elements_start_index, hid_t const datatype_id ) {
#ifndef PERFORMANCE
   // Check whether a group is already opened where the dataset could be appended
   if( active_group_.empty() ) {
      throw std::runtime_error( "Error reserving dataspace. No group is active!" );
   }
   if( datasets_.find( dataspace_name ) != datasets_.end() ) {
      throw std::runtime_error( "Error reserving dataspace. The dataspace already exists!" );
   }
#endif
   // Declare the dataset info struct to store all required information
   Hdf5Dataset dataset;

   /** Link the current active group to the dataset */
   dataset.group_name_ = active_group_;

   // Defines the dataset space global dimensions (copied here to ensure that all data handling is solely done inside of the hdf5_writer class)
   dataset.total_dimensions_ = dataspace_total_dimensions;
   // Defines the local dimensions that are written to the global dataset dataspace
   dataset.local_dimensions_ = dataset_local_dimensions;
   // Defines the start indices for the local dimensions in the global dimensions matrix ( only the cell position is non-zero)
   dataset.start_indices_         = std::vector<hsize_t>( dataset_local_dimensions.size(), 0 );
   dataset.start_indices_.front() = local_elements_start_index;

   // Specific creation type for all datasets that are written with the created dataspace
   dataset.properties_create_ = H5Pcreate( H5P_DATASET_CREATE );

   // Store the dataset data type
   dataset.datatype_ = datatype_id;
   // Create the dataspace for the full simulation data of the given dimensions (NULL marks that the dataspaced is fixed to the given dimensions)
   dataset.dataspace_id_ = H5Screate_simple( dataset.total_dimensions_.size(), dataset.total_dimensions_.data(), NULL );

   /** Add the new dataset to the map */
   datasets_[dataspace_name] = dataset;
}

/**
 * @brief Close all or single dataset/dataspace that was allocated before.
 * @param name Name of the dataset/dataspace that should be closed.
 */
void Hdf5Manager::CloseDataset( std::string const& name ) {
   // If the input string is empty, close all datasets
   if( name.empty() ) {
      // No range-based loop possible here, since we are deleting elements during the loop
      for( auto it = datasets_.begin(); it != datasets_.end(); ) {
         ( *it ).second.Close();
         it = datasets_.erase( it );
      }
   } else {
#ifndef PERFORMANCE
      // Check if the dataset was opened before
      if( datasets_.find( name ) == datasets_.end() ) {
         throw std::runtime_error( "Before closing a dataset/dataspace it must be opened/reserved!" );
      }
#endif
      datasets_[name].Close();
      datasets_.erase( name );
   }
}

/**
 * @brief Close all or single group that have been opened before.
 * @param name Name of the group that should be closed.
 */
void Hdf5Manager::CloseGroup( std::string const& group_name ) {
   // If the input string is empty, close all datasets
   if( group_name.empty() ) {
      // Check if any dataset is still opened, if yes close all
      if( !datasets_.empty() ) {
         CloseDataset();
      }
      // No range-based loop possible here, since we are deleting elements during the loop
      for( auto it = groups_.begin(); it != groups_.end(); ) {
         // erase group from map (calls destructor that closes everything)
         ( *it ).second.Close();
         it = groups_.erase( it );
      }
   } else {
#ifndef PERFORMANCE
      // Check if the group was opened before
      if( groups_.find( group_name ) == groups_.end() ) {
         throw std::runtime_error( "Before closing a group it must be opened!" );
      }
#endif
      // Check if any dataset/dataspace has been opened in that group and close them
      for( auto it = datasets_.begin(); it != datasets_.end(); ) {
         if( ( *it ).second.group_name_ == group_name ) {
            ( *it ).second.Close();
            it = datasets_.erase( it );
         }
      }
      // Then close the desired group and erase it from the map
      groups_[group_name].Close();
      groups_.erase( group_name );
   }
}

/**
 * @brief Close a file which has been opened before.
 */
void Hdf5Manager::CloseFile() {
#ifndef PERFORMANCE
   if( !file_.is_open_ ) {
      throw std::runtime_error( "Cannot close hdf5 file. No file is opened yet.!" );
   }
#endif
   // Check if any group is still opened, if yes close them
   if( !groups_.empty() ) {
      CloseGroup();
   }
   // Close the file
   file_.Close();
}

/**
 * @brief This function is used to get the Hdf5Manager. If no Hdf5Manager exists yet it is created, otherwise the existing writer is passed back. "Singleton Constructor"
 * @return The hdf5_writer instance.
 */
Hdf5Manager& Hdf5Manager::Instance() {
   static Hdf5Manager instance_;
   return instance_;
}