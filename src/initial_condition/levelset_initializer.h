//===--------------------- levelset_initializer.h -------------------------===//
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
#ifndef LEVELSET_INITIALIZER_H
#define LEVELSET_INITIALIZER_H

#include "topology/node.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief The LevelsetInitializer class defines an interface for different
 * initialization methods for the levelset field. It contains the loop structure
 * and the derived classes implement the levelset computation.
 * @note Each bounding box must fully enclose the negative field inside (e.g.,
 * if a plane is provided in 3D, the bounding box must extent to the full
 * negative region and not only the interface).
 */
class LevelsetInitializer {
  // Member of this class
  std::vector<std::array<double, 6>> const bounding_boxes_;
  std::vector<MaterialName> const material_names_;
  double const dimensionalized_node_size_on_level_zero_;
  unsigned int const maximum_level_;

  // Functions for this class only
  bool IsPointInsideBoundingBox(std::array<double, 3> const &point) const;

protected:
  // Protected member
  std::vector<std::string> const spatial_variable_names_ = {"x", "y", "z"};

  // Functions that need to be implemented by the derived classes
  virtual double
  ComputeSignedLevelsetValue(std::array<double, 3> const &point) = 0;
  virtual std::string GetTypeLogData(unsigned int const indent) const = 0;

  // protected default constructor (can only be called from derived classes)
  explicit LevelsetInitializer(
      std::vector<std::array<double, 6>> const bounding_boxes,
      std::vector<MaterialName> const &material_names,
      double const node_size_on_level_zero, unsigned int const maximum_level);

public:
  virtual ~LevelsetInitializer() = default;
  LevelsetInitializer(LevelsetInitializer const &) = delete;
  LevelsetInitializer &operator=(LevelsetInitializer const &) = delete;
  LevelsetInitializer(LevelsetInitializer &&) = delete;
  LevelsetInitializer &operator=(LevelsetInitializer &&) = delete;

  // Public function that can be called from outside
  void GetInitialLevelset(
      nid_t const node_id,
      double (&initial_levelset)[CC::TCX()][CC::TCY()][CC::TCZ()]);
  std::vector<MaterialName> GetInitialMaterials(nid_t const node_id);
  std::string GetLogData(unsigned int const indent) const;
};

#endif // LEVELSET_INITIALIZER_H
