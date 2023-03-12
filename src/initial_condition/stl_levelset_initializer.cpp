//===------------------ stl_levelset_initializer.cpp ----------------------===//
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
#include "stl_levelset_initializer.h"

/**
 * @brief Constructs a stl levelset initializer with levelset initialization
 * parameters given as input.
 * @param dimensional_levelset_initializer_data Map containing all data for the
 * levelset initializer.
 * @param bounding_boxes The bounding boxes inside which an interface can be
 * found.
 * @param material_names Names of the materials.
 * @param node_size_on_level_zero Size of node on level zero.
 * @param maximum_level Maximum level of the simulation.
 */
StlLevelsetInitializer::StlLevelsetInitializer(
    std::string const &stl_filename,
    std::vector<std::array<double, 6>> const &bounding_boxes,
    std::vector<MaterialName> const &material_names,
    double const node_size_on_level_zero, unsigned int const maximum_level)
    : LevelsetInitializer(bounding_boxes, material_names,
                          node_size_on_level_zero, maximum_level),
      stl_filename_(stl_filename),
      stl_triangles_(StlUtilities::ReadStl(stl_filename_)) {
  /* Empty besides initializer list*/
}

/**
 * @brief Computes the level value and appropriate sign for the given point.
 * @param point The point for which the levelset value and sign should be
 * obtained.
 * @return The signed levelset value.
 */
double StlLevelsetInitializer::ComputeSignedLevelsetValue(
    std::array<double, 3> const &point) {
  double levelset = std::numeric_limits<double>::max();
  for (StlUtilities::Triangle const &triangle : stl_triangles_) {
    StlUtilities::Voxelization(triangle, point, levelset);
  }
  return levelset;
}

/**
 * @brief Gives the data for logging for this appropriate class.
 * @param indent Number of white spaces used at the beginning of each line for
 * the logging information.
 * @return string with logging information.
 */
std::string
StlLevelsetInitializer::GetTypeLogData(unsigned int const indent) const {
  // string that is returned
  std::string log_string;
  // Name of the levelset initializer
  log_string += StringOperations::Indent(indent) + "Type               : STL\n";
  // Function
  log_string += StringOperations::Indent(indent) +
                "Filename           : " + stl_filename_ + "\n";
  // Number of triangles
  log_string += StringOperations::Indent(indent) + "Number of triangles: " +
                std::to_string(stl_triangles_.size()) + "\n";

  return log_string;
}
