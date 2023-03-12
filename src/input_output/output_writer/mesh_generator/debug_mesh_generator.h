//===---------------------- debug_mesh_generator.h ------------------------===//
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
#ifndef DEBUG_MESH_GENERATOR_H
#define DEBUG_MESH_GENERATOR_H

#include "input_output/output_writer/mesh_generator.h"
#include <hdf5.h>

/**
 * @brief The DebugMeshGenerator generates a mesh for the output (currently xdmf
 * + hdf5) for all nodes on all levels. It provides the functionality of a
 * non-ambiguous mesh printing the full topology information. This includes all
 * nodes on all levels including the halo cells. Non-ambiguity is provided by
 * placing a gap between different blocks and levels.
 */
class DebugMeshGenerator : public MeshGenerator {

  // Variable specification from the base class
  using MeshGenerator::dimensionalized_node_size_on_level_zero_;
  using MeshGenerator::topology_;
  using MeshGenerator::tree_;

  // Additional parameter for the offset of different levels in the output
  double const z_coordinates_offset_;

  // virtual functions required from the base class to compute data to hdf5 file
  void DoComputeVertexIDs(
      std::vector<unsigned long long int> &vertex_ids) const override;
  void DoComputeVertexCoordinates(
      std::vector<double> &vertex_coordinates) const override;

  // virtual dimension functions required from base class
  std::vector<std::reference_wrapper<Node const>>
  DoGetLocalNodes() const override;
  hsize_t DoGetGlobalNumberOfCells() const override;
  hsize_t DoGetLocalNumberOfCells() const override;
  hsize_t DoGetLocalCellsStartIndex() const override;
  std::vector<hsize_t>
  DoGetGlobalDimensionsOfVertexCoordinates() const override;
  std::vector<hsize_t> DoGetLocalDimensionsOfVertexCoordinates() const override;
  hsize_t DoGetLocalVertexCoordinatesStartIndex() const override;

public:
  DebugMeshGenerator() = delete;
  explicit DebugMeshGenerator(
      TopologyManager const &topology, Tree const &flower,
      double const dimensionalized_node_size_on_level_zero,
      unsigned int const number_of_z_nodes_on_level_zero);
  virtual ~DebugMeshGenerator() = default;
  DebugMeshGenerator(DebugMeshGenerator const &) = delete;
  DebugMeshGenerator &operator=(DebugMeshGenerator const &) = delete;
  DebugMeshGenerator(DebugMeshGenerator &&) = delete;
  DebugMeshGenerator &operator=(DebugMeshGenerator &&) = delete;
};

#endif // DEBUG_MESH_GENERATOR_H
