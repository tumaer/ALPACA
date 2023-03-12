//===------------------- stl_levelset_initializer.h -----------------------===//
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
#ifndef STL_LEVELSET_INITIALIZER_H
#define STL_LEVELSET_INITIALIZER_H

#include "initial_condition/levelset_initializer.h"
#include "input_output/utilities/stl_utilities.h"

/**
 * @brief The StlLevelsetInitializer class allows for a levelset initialization
 * based on a provided STL file. It contains the levelset computation, the loop
 * structure is in the base class.
 * @note The number and size of the triangles, the surface is composed of, must
 * be appropriate for the used mesh size. In case smaller cell sizes are used,
 *       the number of triangles must be increased and the size reduced. The
 * inversely holds for coarser meshes.
 */
class StlLevelsetInitializer : public LevelsetInitializer {
  // Member variables of this class only
  std::string const stl_filename_;

  std::vector<StlUtilities::Triangle> stl_triangles_;

  // Functions required from base class
  double
  ComputeSignedLevelsetValue(std::array<double, 3> const &point) override;
  std::string GetTypeLogData(unsigned int const indent) const override;

public:
  StlLevelsetInitializer() = delete;
  explicit StlLevelsetInitializer(
      std::string const &stl_filename,
      std::vector<std::array<double, 6>> const &bounding_boxes,
      std::vector<MaterialName> const &material_names,
      double const node_size_on_level_zero, unsigned int const maximum_level);
  virtual ~StlLevelsetInitializer() = default;
  StlLevelsetInitializer(StlLevelsetInitializer const &) = delete;
  StlLevelsetInitializer &operator=(StlLevelsetInitializer const &) = delete;
  StlLevelsetInitializer(StlLevelsetInitializer &&) = delete;
  StlLevelsetInitializer &operator=(StlLevelsetInitializer &&) = delete;
};

#endif // STL_LEVELSET_INITIALIZER_H
