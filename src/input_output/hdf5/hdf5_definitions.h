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
#ifndef HDF5_DEFINITIONS_H
#define HDF5_DEFINITIONS_H

#include <hdf5.h>
#include <vector>
#include <string>

/**
 * @brief Enum class for the different access options of the hdf5 file.
 */
enum class Hdf5Access { Read,
                        Write };

/**
 * @brief Struct that provides all information required for accessing (reading/writing) data from/to the hdf5 file. It provides all information that is required
 *        to access the data through datasets and dataspaces. This struct is used in the Hdf5Manager class. Depending on the chosen order of reading/writing
 *        procedures, the struct can only be filled partially.
 *        Naming:
 *          Dataspace: Specifies the global extent of the data (covers all nodes/cells on all ranks).
 *          Dataset: Specifies the set that is linked into the dataspace with correct name and group association (no size definition).
 *          Local Memory Space: The allocated memory space that is written/read locally (on current rank) during one reading/writing process.
 *          Local Hypeslab: Specifies the positions where the local memory space is written/read to/from in the dataset.
 *
 *          Total Dimensions: Specifies the dimensions of the dataspace that is written totally from all ranks (=> Dataspace).
 *          Local Dimensions: Specifies the dimensions each rank writes during one writing/reading process (=> Local Memory Space).
 *          Start Indices: Specifies the positions where the data is written/read to/from (=> local hyperslab).
 */
struct Hdf5Dataset {
   // String holding the group name, where the data is placed into
   std::string group_name_;
   // Dimension definitons for a dataset/dataspace (all need the same rank)
   std::vector<hsize_t> total_dimensions_;
   std::vector<hsize_t> local_dimensions_;
   std::vector<hsize_t> start_indices_;
   std::vector<hsize_t> chunk_;
   std::vector<hsize_t> count_;
   // Identifier for opened and/or created datasets
   hid_t datatype_           = -1;
   hid_t properties_create_  = -1;
   hid_t dataspace_id_       = -1;
   hid_t dataset_id_         = -1;
   hid_t local_memory_space_ = -1;
   hid_t local_hyperslab_    = -1;

   void Close() {
      if( properties_create_ != -1 ) H5Pclose( properties_create_ );
      if( local_hyperslab_ != -1 ) H5Sclose( local_hyperslab_ );
      if( local_memory_space_ != -1 ) H5Sclose( local_memory_space_ );
      if( dataset_id_ != -1 ) H5Dclose( dataset_id_ );
      if( dataspace_id_ != -1 ) H5Sclose( dataspace_id_ );
   }
};

/**
 * @brief Struct holding the information that is required to open/create groups in the hdf5 file.
 */
struct Hdf5Group {
   hid_t id_         = -1;
   hid_t properties_ = -1;

   void Close() {
      if( properties_ != -1 ) H5Pclose( properties_ );
      if( id_ != -1 ) H5Gclose( id_ );
   }
};

/**
 * @brief Struct holding the information that is required to open/create hdf5 files.
 */
struct Hdf5File {
   bool is_open_ = false;
   Hdf5Access access_type_;
   hid_t id_;
   hid_t properties_;

   void Close() {
      if( properties_ != -1 ) H5Pclose( properties_ );
      if( id_ != -1 ) {
         H5Fclose( id_ );
         is_open_ = false;
      }
   }
};

#endif// HDF5_DEFINITIONS_H