//===---------------------- mesh_generator.cpp ----------------------------===//
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
#include "mesh_generator.h"
#include "input_output/utilities/xdmf_utilities.h"

/**
 * @brief Constructor for a generic mesh generator to be called from derived
 * classes.
 * @param topology_manager Manager storing all topology information of all
 * ranks.
 * @param flower The tree for the node informations on the current rank.
 * @param node_size_on_level_zero Already dimensionalized size of a single node
 * on level zero.
 */
MeshGenerator::MeshGenerator(
    TopologyManager const &topology_manager, Tree const &flower,
    double const dimensionalized_node_size_on_level_zero)
    : topology_(topology_manager), tree_(flower),
      dimensionalized_node_size_on_level_zero_(
          dimensionalized_node_size_on_level_zero) {
  /** Empty constructor besides initializer list */
}

/**
 * @brief Returns the string for the topology (vertex ids of cells).
 * @param filename Name of the .hdf5 file where the actual data is found.
 * @param group_name Name of the group the topology is written into.
 * @return Attribute string for the topology.
 */
std::string
MeshGenerator::GetXdmfTopologyString(std::string const &filename,
                                     std::string const &group_name) const {
  hsize_t const global_number_cells = GetGlobalNumberOfCells();
  std::string const data_item(
      "<DataItem NumberType=\"Int\" Format=\"HDF\" Dimensions=\"" +
      std::to_string(global_number_cells) + " 8\"> " + filename + ":/" +
      group_name + "/" + vertex_ids_name_ + " </DataItem>\n");
  return XdmfUtilities::TopologyString(data_item, global_number_cells);
}

/**
 * @brief Returns the string used for the geometry (vertex coordinates).
 * @param filename Name of the .hdf5 file where the actual data is found.
 * @param group_name Name of the group the geometry is written into.
 * @return Attribute string for the geometry.
 */
std::string
MeshGenerator::GetXdmfGeometryString(std::string const &filename,
                                     std::string const &group_name) const {
  hsize_t const global_number_vertices =
      GetGlobalDimensionsOfVertexCoordinates().front();
  std::string const data_item("<DataItem Format=\"HDF\" NumberType=\"Float\" "
                              "Precision=\"8\" Dimensions=\"" +
                              std::to_string(global_number_vertices) +
                              " 3\"> " + filename + ":/" + group_name + "/" +
                              vertex_coordinates_name_ + " </DataItem>\n");
  return XdmfUtilities::GeometryString(data_item, global_number_vertices);
}

/**
 * @brief Appends the vertex IDs for the specific mesh generator to the hdf5
 * file. The vertex IDs represent the cells, where one cell is built by 8
 * vertices.
 * @param vertex_ids Vector where all vertex IDs are written into (indirect
 * return).
 */
void MeshGenerator::ComputeVertexIDs(
    std::vector<unsigned long long int> &vertex_ids) const {
  // Compute the vertex IDs in the derived class
  DoComputeVertexIDs(vertex_ids);
}

/**
 * @brief Appends the vertex coordinates for the specific mesh generator to the
 * hdf5 file.
 * @param vertex_coordinates Vector where all vertex coordinates are written
 * into (indirect return).
 */
void MeshGenerator::ComputeVertexCoordinates(
    std::vector<double> &vertex_coordinates) const {
  // Compute the vertex coordinates in the derived class
  DoComputeVertexCoordinates(vertex_coordinates);
}
