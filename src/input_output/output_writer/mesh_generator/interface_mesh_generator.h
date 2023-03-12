//===------------------- interface_mesh_generator.h -----------------------===//
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
#ifndef INTERFACE_MESH_GENERATOR_H
#define INTERFACE_MESH_GENERATOR_H

#include "input_output/output_writer/mesh_generator.h"
#include <hdf5.h>

/**
 * @brief The InterfaceMeshGenerator generates a mesh for the output (currently
 * xdmf + hdf5) using only the interface leaves. It provides the functionality
 * of a (currently) ambiguous mesh (vertex IDs and coordinates) for all nodes
 * containing an interface. Ambiguity is generated, since in parallel writing,
 * different ranks generate IDs and coordinates for the same vertex (at block
 * borders). For remove the ambiguity a vertex filter should be implemented or a
 * complete different routine must be provided. The final mesh represents all
 * nodes containing interfaces (always finest level). Only leaf nodes are
 * written.
 */
class InterfaceMeshGenerator : public MeshGenerator {

  // Variable specification from the base class
  using MeshGenerator::dimensionalized_node_size_on_level_zero_;
  using MeshGenerator::topology_;
  using MeshGenerator::tree_;

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
  InterfaceMeshGenerator() = delete;
  explicit InterfaceMeshGenerator(
      TopologyManager const &topology, Tree const &flower,
      double const dimensionalized_node_size_on_level_zero);
  virtual ~InterfaceMeshGenerator() = default;
  InterfaceMeshGenerator(InterfaceMeshGenerator const &) = delete;
  InterfaceMeshGenerator &operator=(InterfaceMeshGenerator const &) = delete;
  InterfaceMeshGenerator(InterfaceMeshGenerator &&) = delete;
  InterfaceMeshGenerator &operator=(InterfaceMeshGenerator &&) = delete;
};

#endif // INTERFACE_MESH_GENERATOR_H
