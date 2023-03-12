//===----------------- standard_mpi_mesh_generator.h ----------------------===//
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
#ifndef STANDARD_MPI_MESH_GENERATOR_H
#define STANDARD_MPI_MESH_GENERATOR_H

#include "input_output/output_writer/mesh_generator.h"
#include <hdf5.h>

/**
 * @brief The StandardMpiMeshGenerator generates a mesh for the output
 * (currently xdmf + hdf5) using a mpi routine.
 *
 *        It provides the functionality of a non-ambiguous mesh (vertex IDs and
 * coordinates) using a mpi-routine. If the mpi-routine is not used double
 * placed vertices will be present in the output. The mesh represents the
 * current multi-resolution situation including jumps between different blocks.
 * Only leaf nodes are written.
 */
class StandardMpiMeshGenerator : public MeshGenerator {

  // Variable specification from the base class
  using MeshGenerator::dimensionalized_node_size_on_level_zero_;
  using MeshGenerator::topology_;
  using MeshGenerator::tree_;

  // self defined member variables
  bool const mpi_filtering_active_;

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

  // local function for vertex filtering using mpi routing
  void FilterVertexIDs(
      std::vector<unsigned long long int> &vertex_ids,
      std::vector<unsigned long long int> const &leave_offset) const;

public:
  StandardMpiMeshGenerator() = delete;
  explicit StandardMpiMeshGenerator(
      TopologyManager const &topology, Tree const &flower,
      double const dimensionalized_node_size_on_level_zero,
      bool const mpi_filtering_active);
  virtual ~StandardMpiMeshGenerator() = default;
  StandardMpiMeshGenerator(StandardMpiMeshGenerator const &) = delete;
  StandardMpiMeshGenerator &
  operator=(StandardMpiMeshGenerator const &) = delete;
  StandardMpiMeshGenerator(StandardMpiMeshGenerator &&) = delete;
  StandardMpiMeshGenerator &operator=(StandardMpiMeshGenerator &&) = delete;
};

#endif // STANDARD_MPI_MESH_GENERATOR_H
