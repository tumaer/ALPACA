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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#ifndef HDF5_MANAGER_H
#define HDF5_MANAGER_H

#include <hdf5.h>
#include <string>
#include <vector>
#include <stdexcept>
#include <unordered_map>

#include "input_output/hdf5/hdf5_definitions.h"

/**
 * @brief A light-weight hdf5 file writer class to write output to an hdf5 file.
 * @note Singleton.
 */
class Hdf5Manager {

   // Member variables needed for the file to be written to
   Hdf5File file_;
   // Member variables used for a group definition
   std::string active_group_;
   std::unordered_map<std::string, Hdf5Group> groups_;
   // Member variables for specification of a single dataset
   std::unordered_map<std::string, Hdf5Dataset> datasets_;

   // Constructor called from the singleton public constructor
   explicit Hdf5Manager();

public:
   // Singleton constructor
   static Hdf5Manager& Instance();

   ~Hdf5Manager();
   Hdf5Manager( Hdf5Manager const& ) = delete;
   Hdf5Manager& operator=( Hdf5Manager const& ) = delete;
   Hdf5Manager( Hdf5Manager&& )                 = delete;
   Hdf5Manager& operator=( Hdf5Manager&& ) = delete;

   // Functions to open, close files, groups, datasets and dataspaces
   void OpenFile( std::string const& filename, Hdf5Access const access_type = Hdf5Access::Write );
   void CloseFile();
   void OpenGroup( std::string const& group_name );
   void ActivateGroup( std::string const& name );
   void CloseGroup( std::string const& group_name = "" );
   void OpenDatasetForWriting( std::string const& dataset_name, std::vector<hsize_t> const& dataspace_total_dimensions, std::vector<hsize_t> const& dataset_local_dimensions, hsize_t const local_elements_start_index, hid_t const datatype_id = H5T_NATIVE_DOUBLE );
   void OpenDatasetForReading( std::string const& dataset_name, std::vector<hsize_t> const& dataset_local_dimensions, hid_t const datatype_id = H5T_NATIVE_DOUBLE );
   void ReserveDataspace( std::string const& dataspace_name, std::vector<hsize_t> const& dataspace_total_dimensions, std::vector<hsize_t> const& dataset_local_dimensions, hsize_t const local_elements_start_index, hid_t const datatype_id );
   void CloseDataset( std::string const& dataset_name = "" );
   // Gives the extent of a dataset
   hsize_t GetDatasetExtent( std::string const& dataset_name ) const;

   /**
    * @brief Writes a single scalar attribute to the currently opened group.
    * @param scalar_name Name of the scalar that is written.
    * @param value Value of the scalar that is written into the file.
    * @param datatype_id Hdf5 identifier which datatype should be used for the dataset (e.g., H5T_NATIVE_ULLONG, H5T_NATIVE_DOUBLE).
    *
    * @tparam Type of values to be written.
    * @note Attributes do not represent actual datasets. They are used to provide meta data for the group/dataset.
    */
   template<typename T>
   void WriteAttributeScalar( std::string const& scalar_name, T const value, hid_t const datatype_id = H5T_NATIVE_DOUBLE ) const {
#ifndef PERFORMANCE
      if( active_group_.empty() ) {
         throw std::runtime_error( "Before write a single scalar to a group you must activate one!" );
      }
#endif
      // Get the group data
      Hdf5Group const& group = groups_.at( active_group_ );
      // Write to the group
      hid_t data_type = H5Screate( H5S_SCALAR );
      hid_t data      = H5Acreate2( group.id_, scalar_name.c_str(), datatype_id, data_type, H5P_DEFAULT, H5P_DEFAULT );
      H5Awrite( data, datatype_id, &value );
      H5Aclose( data );
      H5Sclose( data_type );
   }

   /**
    * @brief Writes data into an already opened dataset. It is part of a two-phase writing procedure:
    *        1. OpenDatasetForWriting (allocates the memory for the dataset to be written)
    *        2. WriteDataset (selects the hyperslab position and writes data into the
    *        Therefore, before calling this function, the dataset must always be opened.
    * @param dataset_name Name of the dataset that is written (must conincide with the name used to open).
    * @param buffer Pointer to the CONTIGUOUS buffer that is written.
    *
    * @tparam BufferType Buffer type that is written.
    */
   template<typename BufferType>
   void WriteDataset( std::string const& dataset_name, BufferType const* buffer ) {
#ifndef PERFORMANCE
      // Check if the dataset was opened before
      if( datasets_.find( dataset_name ) == datasets_.end() ) {
         throw std::logic_error( "Before writing to a dataset it must be opened!" );
      }
#endif
      // Get the correct dataset info
      Hdf5Dataset& dataset   = datasets_[dataset_name];
      Hdf5Group const& group = groups_[dataset.group_name_];
      // Select the correct hyperslab
      H5Sselect_hyperslab( dataset.local_hyperslab_, H5S_SELECT_SET, dataset.start_indices_.data(), NULL, dataset.count_.data(), dataset.local_dimensions_.data() );
      // Write the local dataset to the given group
      H5Dwrite( dataset.dataset_id_, dataset.datatype_, dataset.local_memory_space_, dataset.local_hyperslab_, group.properties_, buffer );
      // Increment the dataset start index for the next writing process
      dataset.start_indices_.front()++;
   }

   /**
    * @brief Writes data into an already allocated dataspace. It is part of a two-phase writing procedure:
    *           1. ReserveDataspace (allocates the memory for the dataspace to be written)
    *           2. WriteDatasetToDataspace (creates the dataset, selects the hyperslab position and writes data to the dataset)
    *        Therefore, before calling this function, the dataspace must always be reserved. During the function call the dataset is opened and written.
    * @param dataspace_name Name of the dataspace used to write the dataset (must coincide with name that is used to reserve dataspace).
    * @param dataset_name Name of the dataset that is written.
    * @param buffer Pointer to the CONTIGUOUS buffer that is written.
    *
    * @tparam BufferType Buffer type that is written.
    */
   template<typename BufferType>
   void WriteDatasetToDataspace( std::string const& dataspace_name, std::string const& dataset_name, BufferType const* buffer ) const {
#ifndef PERFORMANCE
      if( datasets_.find( dataspace_name ) == datasets_.end() ) {
         throw std::runtime_error( "Before writing a datset to dataspace, a dataspace must be reserved!" );
      }
#endif
      // Get the correct group information
      Hdf5Dataset const& dataset = datasets_.at( dataspace_name );
      Hdf5Group const& group     = groups_.at( dataset.group_name_ );

      // Creates a new dataset and links it into the file/group
      hid_t const dataset_id = H5Dcreate2( group.id_, dataset_name.c_str(), dataset.datatype_, dataset.dataspace_id_, H5P_DEFAULT, dataset.properties_create_, H5P_DEFAULT );
      // Create the local memory and hyperslab space required for writing the data (NULL marks that dataset cannot grow until it is closed)
      hid_t const local_memory_space = H5Screate_simple( dataset.local_dimensions_.size(), dataset.local_dimensions_.data(), NULL );
      hid_t const local_hyperslab    = H5Dget_space( dataset_id );
      H5Sselect_hyperslab( local_hyperslab, H5S_SELECT_SET, dataset.start_indices_.data(), NULL, dataset.local_dimensions_.data(), NULL );
      // Write the local dataset to the given group
      H5Dwrite( dataset_id, dataset.datatype_, local_memory_space, local_hyperslab, group.properties_, buffer );
      // close all local reserved data
      H5Sclose( local_hyperslab );
      H5Sclose( local_memory_space );
      H5Dclose( dataset_id );
   }

   /**
    * @brief Reads a single attribute scalar value from the currently opened group.
    * @param scalar_name Name of the scalar that is read.
    * @param value Value where the read scalar data is read into.
    * @param datatype_id Hdf5 identifier which datatype should be used for the dataset (e.g., H5T_NATIVE_ULLONG, H5T_NATIVE_DOUBLE).
    *
    * @tparam Type of values to be written.
    */
   template<typename T>
   void ReadAttributeScalar( std::string const& scalar_name, T& value, hid_t const datatype_id = H5T_NATIVE_DOUBLE ) const {
#ifndef PERFORMANCE
      if( active_group_.empty() ) {
         throw std::runtime_error( "Before read a single scalar from a group you must activate one!" );
      }
#endif
      // Get the group data
      Hdf5Group const& group = groups_.at( active_group_ );
      // Read from the group
      hid_t scalar_attribute = H5Aopen( group.id_, scalar_name.c_str(), H5P_DEFAULT );
      H5Aread( scalar_attribute, datatype_id, &value );
      H5Aclose( scalar_attribute );
   }

   /**
    * @brief Reads a full (contiguous) dataset in one chunk.
    * @param dataset_name Name of the dataset that should be read.
    * @param buffer Buffer where the read dataset is stored into.
    * @param datatype_id Hdf5 identifier which datatype should be used for the dataset (e.g., H5T_NATIVE_ULLONG, H5T_NATIVE_DOUBLE).
    *
    * @tparam BufferType Type of the buffer where data is stored.
    * @note No sanity check is done that the provided buffer is large enough to store the data. Use GetDatasetExtent function to check the (contiguous) size of dataset.
    */
   template<typename BufferType>
   void ReadFullDataset( std::string const& dataset_name, BufferType* buffer, hid_t const datatype_id = H5T_NATIVE_DOUBLE ) const {
#ifndef PERFORMANCE
      if( active_group_.empty() ) {
         throw std::runtime_error( "Before reading a full dataset from a hdf5 file, a group must be opened and activated!" );
      }
#endif
      // Get the correct group info
      Hdf5Group const& group = groups_.at( active_group_ );

      /** Open the dataset for reading and select the entire local hyperslab */
      hid_t dataset_id      = H5Dopen2( group.id_, dataset_name.c_str(), H5P_DEFAULT );
      hid_t local_hyperslab = H5Dget_space( dataset_id );
      // read the data from file and store it in the contiguous buffer
      H5Dread( dataset_id, datatype_id, H5S_ALL, H5S_ALL, group.properties_, buffer );
      // close the hyperslab and dataset
      H5Sclose( local_hyperslab );
      H5Dclose( dataset_id );
   }

   /**
    * @brief Reads a single dataset with a given offset from the opened dataset. It is part of a two-phase reading procedure:
    *           1. OpenDatasetForReading (allocates the memory for the dataset to be read)
    *           2. ReadDataset (selects the hyperslab position and read data from the dataset)
    *        Therefore, always open the dataset for reading before.
    * @param dataset_name Name of the dataset that should be read.
    * @param buffer Buffer where the read dataset is stored into.
    * @param dataset_offset Offset position where the buffer should start the reading procedure.
    * @tparam BufferType Type of the buffer where data is stored.
    */
   template<typename BufferType>
   void ReadDataset( std::string const& dataset_name, BufferType* buffer, hsize_t const dataset_offset ) {
#ifndef PERFORMANCE
      // Check if the dataset was opened before
      if( datasets_.find( dataset_name ) == datasets_.end() ) {
         throw std::logic_error( "Before writing to a dataset it must be opened!" );
      }
#endif
      // Get the correct dataset info
      Hdf5Dataset& dataset   = datasets_[dataset_name];
      Hdf5Group const& group = groups_[dataset.group_name_];
      // Change the offset to the desired enty
      dataset.start_indices_.front() = dataset_offset;
      H5Sselect_hyperslab( dataset.local_hyperslab_, H5S_SELECT_SET, dataset.start_indices_.data(), NULL, dataset.count_.data(), dataset.local_dimensions_.data() );
      H5Dread( dataset.dataset_id_, dataset.datatype_, dataset.local_memory_space_, dataset.local_hyperslab_, group.properties_, buffer );
   }
};

#endif// HDF5_MANAGER_H
