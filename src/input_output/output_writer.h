//===------------------------ output_writer.h -----------------------------===//
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
#ifndef OUTPUT_WRITER_H
#define OUTPUT_WRITER_H

#include <map>

#include "input_output/hdf5/hdf5_manager.h"
#include "input_output/input_reader/multi_resolution_reader/multi_resolution_reader.h"
#include "input_output/input_reader/output_reader/output_reader.h"
#include "input_output/output_writer/mesh_generator.h"
#include "input_output/output_writer/output_quantity.h"
#include "materials/material_manager.h"

/**
 * @brief The OutputWriter class handles the output to the filesystem in XDMF +
 * HDF5 file format for ParaView. OutputWriter must not change any data. It
 * provides the functionality that a xdmf and hdf5 file is written. Therefore,
 * the mesh is generated in the specific format depending on the desired output
 * type (standard, interface, debug). Furthermore, an xdmf file is written that
 * provides the direct access to all files for the simulation in a time series.
 *        The xdmf file can be used with the Xdmf2 and Xdmf3 reader in ParaView.
 */
class OutputWriter {
  // Specification that this variable is used from the base class to avoid
  // shadowing
  unsigned int const number_of_materials_;

  // Instance of the hdf5 file writer (cannot be const due to variable changes
  // during simulation, singleton allows constness of OutputWriter)
  Hdf5Manager &hdf5_manager_;

  // The different mesh generators for the standard, debug and interface output
  // (vertexIDs and coordinates generation)
  std::unique_ptr<MeshGenerator const> const standard_mesh_generator_;
  std::unique_ptr<MeshGenerator const> const debug_mesh_generator_;
  std::unique_ptr<MeshGenerator const> const interface_mesh_generator_;

  // vector containing all output quantities used for the output (unique_ptr
  // since it is the base class) In general both vectors can be used in one
  // single, but in the writing process a distinction is made for the debug
  // output. Therefore, separated.
  std::vector<std::unique_ptr<OutputQuantity const>> const
      material_output_quantities_;
  std::vector<std::unique_ptr<OutputQuantity const>> const
      interface_output_quantities_;

  // Map that provides the dimensionalization information of all quantities
  // (allows single vector allocations for quantities with same dimensions)
  std::map<std::array<unsigned int, 2>, std::vector<unsigned int>> const
      material_quantities_dimension_map_;
  std::map<std::array<unsigned int, 2>, std::vector<unsigned int>> const
      interface_quantities_dimension_map_;

  // local functions to write the hdf5 and xdmf files
  void WriteHdf5File(double const output_time, std::string const &hdf5_filename,
                     MeshGenerator const &mesh_generator,
                     OutputType const output_type) const;
  void WriteXdmfTimeStepFile(double const output_time,
                             std::string const &hdf5_filename,
                             MeshGenerator const &mesh_generator,
                             OutputType const output_type) const;
  void AppendToXdmfTimeSeriesFile(
      double const output_time, std::string const &hdf5_filename,
      MeshGenerator const &mesh_generator,
      std::string const &time_series_filename_without_extension,
      OutputType const output_type) const;
  std::string XdmfSpatialDataInformation(double const output_time,
                                         std::string const &hdf5_short_filename,
                                         MeshGenerator const &mesh_generator,
                                         OutputType const output_type) const;

  // local factory functions
  std::map<std::array<unsigned int, 2>, std::vector<unsigned int>>
  ComputeDimensionMap(std::vector<std::unique_ptr<OutputQuantity const>> const
                          &quantities) const;

public:
  OutputWriter() = delete;
  explicit OutputWriter(
      std::unique_ptr<MeshGenerator const> standard_mesh_generator,
      std::unique_ptr<MeshGenerator const> debug_mesh_generator,
      std::unique_ptr<MeshGenerator const> interface_mesh_generator,
      std::vector<std::unique_ptr<OutputQuantity const>>
          material_output_quantities,
      std::vector<std::unique_ptr<OutputQuantity const>>
          interface_output_quantities,
      unsigned int const number_of_materials);
  ~OutputWriter() = default;
  OutputWriter(OutputWriter const &) = delete;
  OutputWriter &operator=(OutputWriter const &) = delete;
  OutputWriter(OutputWriter &&) = delete;
  OutputWriter &operator=(OutputWriter &&) = delete;

  // Function to write the actual output files
  void WriteOutput(
      OutputType const output_type, double const output_time,
      std::string const &filename_without_extension,
      std::string const &time_series_filename_without_extension = "") const;
  // Function to initialize and finalize the time series files
  void InitializeTimeSeriesFile(
      std::string const &time_series_filename_without_extension) const;
  void FinalizeTimeSeriesFile(
      std::string const &time_series_filename_without_extension) const;
};

#endif // OUTPUT_WRITER_H
