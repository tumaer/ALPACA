//===----------------------- hdf5_definitions.h ---------------------------===//
//
//                                 ALPACA
//
// Part of ALPACA, under the GNU General Public License as published by
// the Free Software Foundation version 3.
// SPDX-License-Identifier: GPL-3.0-only
//
// If using this code in an academic setting, please cite the following:
// @article{hoppe2022parallel,
//  title={A parallel modular computing environment for three-dimensional
//  multiresolution simulations of compressible flows},
//  author={Hoppe, Nils and Adami, Stefan and Adams, Nikolaus A},
//  journal={Computer Methods in Applied Mechanics and Engineering},
//  volume={391},
//  pages={114486},
//  year={2022},
//  publisher={Elsevier}
// }
//
//===----------------------------------------------------------------------===//
#ifndef HDF5_DEFINITIONS_H
#define HDF5_DEFINITIONS_H

#include <hdf5.h>
#include <string>
#include <vector>

/**
 * @brief Enum class for the different access options of the hdf5 file.
 */
enum class Hdf5Access { Read, Write };

/**
 * @brief Struct that provides all information required for accessing
 * (reading/writing) data from/to the hdf5 file. It provides all information
 * that is required to access the data through datasets and dataspaces. This
 * struct is used in the Hdf5Manager class. Depending on the chosen order of
 * reading/writing procedures, the struct can only be filled partially. Naming:
 *          Dataspace: Specifies the global extent of the data (covers all
 * nodes/cells on all ranks). Dataset: Specifies the set that is linked into the
 * dataspace with correct name and group association (no size definition). Local
 * Memory Space: The allocated memory space that is written/read locally (on
 * current rank) during one reading/writing process. Local Hypeslab: Specifies
 * the positions where the local memory space is written/read to/from in the
 * dataset.
 *
 *          Total Dimensions: Specifies the dimensions of the dataspace that is
 * written totally from all ranks (=> Dataspace). Local Dimensions: Specifies
 * the dimensions each rank writes during one writing/reading process (=> Local
 * Memory Space). Start Indices: Specifies the positions where the data is
 * written/read to/from (=> local hyperslab).
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
  hid_t datatype_ = -1;
  hid_t properties_create_ = -1;
  hid_t dataspace_id_ = -1;
  hid_t dataset_id_ = -1;
  hid_t local_memory_space_ = -1;
  hid_t local_hyperslab_ = -1;

  void Close() {
    if (properties_create_ != -1)
      H5Pclose(properties_create_);
    if (local_hyperslab_ != -1)
      H5Sclose(local_hyperslab_);
    if (local_memory_space_ != -1)
      H5Sclose(local_memory_space_);
    if (dataset_id_ != -1)
      H5Dclose(dataset_id_);
    if (dataspace_id_ != -1)
      H5Sclose(dataspace_id_);
  }
};

/**
 * @brief Struct holding the information that is required to open/create groups
 * in the hdf5 file.
 */
struct Hdf5Group {
  hid_t id_ = -1;
  hid_t properties_ = -1;

  void Close() {
    if (properties_ != -1)
      H5Pclose(properties_);
    if (id_ != -1)
      H5Gclose(id_);
  }
};

/**
 * @brief Struct holding the information that is required to open/create hdf5
 * files.
 */
struct Hdf5File {
  bool is_open_ = false;
  Hdf5Access access_type_;
  hid_t id_;
  hid_t properties_;

  void Close() {
    if (properties_ != -1)
      H5Pclose(properties_);
    if (id_ != -1) {
      H5Fclose(id_);
      is_open_ = false;
    }
  }
};

#endif // HDF5_DEFINITIONS_H
