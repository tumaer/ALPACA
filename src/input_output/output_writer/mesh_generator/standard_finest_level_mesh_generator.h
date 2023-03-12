//===------------- standard_finest_level_mesh_generator.h -----------------===//
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
#ifndef FINEST_LEVEL_MESH_GENERATOR_H
#define FINEST_LEVEL_MESH_GENERATOR_H

#include "input_output/output_writer/mesh_generator.h"
#include <hdf5.h>

/**
 * @brief The StandardFinestLevelMeshGenerator generates a mesh for the output
 * (currently xdmf + hdf5) using a virtual mapping of the mesh on the finest
 * level.
 *
 *        It provides the functionality of a non-ambiguous mesh (vertex IDs and
 * coordinates) using the currently finest level as a reference. All blocks are
 * virtually mapped onto the current maximum level. Therefore, the final hdf5
 * file is overloaded with vertex coordinates of the finest level, where a lot
 * of are not used in the reading, e.g. using ParaView. Nevertheless, the final
 * mesh represents the current multi-resolution situation including jumps
 * between different blocks. Only leaf nodes are written.
 *
 *        Example for specification of vertex IDs:
 *
 *        Actual mesh situation:                         Virtual mapping on
 * finest level
 *
 *         4____5________6_______________7 8____9____10___11___12___13___14___15
 *         |____|________|_______________| |____|____|____|____|____|____|____|
 *         0    1        2               3                0    1    2    3    4
 * 5    6    7
 *
 *        Final Output:
 *
 *         8____9________11______________15
 *         |____|________|_______________|
 *         0    1        3               7
 *
 *        Coordinates written to hdf5: all from 0 to 15, even only 8 are used,
 * to provide direct mapping from vertex IDs to coordinates (required for
 * reading)
 *
 */
class StandardFinestLevelMeshGenerator : public MeshGenerator {

  // Variable specification from the base class
  using MeshGenerator::dimensionalized_node_size_on_level_zero_;
  using MeshGenerator::topology_;
  using MeshGenerator::tree_;

  // Self defined variables
  std::array<unsigned int, 3> const number_of_nodes_on_level_zero_;

  // local functions required for the computation
  std::array<unsigned long long int, 3>
  GlobalNumberOfVerticesPerDimension() const;
  unsigned long long int TotalNumberOfVerticesOnFinestLevel() const;
  std::array<unsigned long long int, 2>
  GetLocalNumberOfVerticesAndStartIndex() const;

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
  StandardFinestLevelMeshGenerator() = delete;
  explicit StandardFinestLevelMeshGenerator(
      TopologyManager const &topology_manager, Tree const &flower,
      double const dimensionalized_node_size_on_level_zero,
      std::array<unsigned int, 3> const number_of_nodes_on_level_zero);
  virtual ~StandardFinestLevelMeshGenerator() = default;
  StandardFinestLevelMeshGenerator(StandardFinestLevelMeshGenerator const &) =
      delete;
  StandardFinestLevelMeshGenerator &
  operator=(StandardFinestLevelMeshGenerator const &) = delete;
  StandardFinestLevelMeshGenerator(StandardFinestLevelMeshGenerator &&) =
      delete;
  StandardFinestLevelMeshGenerator &
  operator=(StandardFinestLevelMeshGenerator &&) = delete;
};

#endif // FINEST_LEVEL_MESH_GENERATOR_H
