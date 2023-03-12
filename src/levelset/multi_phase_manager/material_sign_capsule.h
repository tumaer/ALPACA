//===---------------------- material_sign_capsule.h -----------------------===//
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
#ifndef MATERIAL_SIGN_CAPSULE_H
#define MATERIAL_SIGN_CAPSULE_H

#include "materials/material_definitions.h"

/**
 * @brief Work-around class to reduce cyclic dependencies and still query the
 * material sign where needed. It is the user's responsibility to call functions
 * only after initialization (currently done in the constructor of the material
 * manager class that holds all information for the two materials).
 */
class MaterialSignCapsule {

private:
  static MaterialName positive_material_;
  static MaterialName negative_material_;

public:
  MaterialSignCapsule() = delete;
  ~MaterialSignCapsule() = default;
  MaterialSignCapsule(MaterialSignCapsule const &) = delete;
  MaterialSignCapsule &operator=(MaterialSignCapsule const &) = delete;
  MaterialSignCapsule(MaterialSignCapsule &&) = delete;
  MaterialSignCapsule &operator=(MaterialSignCapsule &&) = delete;

  MaterialSignCapsule(MaterialName const positive_material,
                      MaterialName const negative_material) {
    positive_material_ = positive_material;
    negative_material_ = negative_material;
  }

  /**
   * @brief Static function to get the positive material material identifier in
   * a single-level-set simulation. $Always MATERIAL ONE in inputfile. This
   * function can be uninitialized if called too early! Must not be called as
   * long as no object is available.$
   * @return Positive material material identifier .
   */
  static inline MaterialName PositiveMaterial() { return positive_material_; }

  /**
   * @brief Static function to get the negative material material identifier in
   * a single-level-set simulation. $Always MATERIAL TWO in inputfile. This
   * function can be uninitialized if called too early! Must not be called as
   * long as no object is available.$
   * @return Negative material material identifier.
   */
  static inline MaterialName NegativeMaterial() { return negative_material_; }

  /**
   * @brief Gives the sign of the given material used in the signed levelset and
   * signed interface tag description.
   * @param material Material of interest.
   * @return return Sign of the material.
   */
  static inline std::int8_t SignOfMaterial(MaterialName const material) {
    return material == PositiveMaterial() ? 1 : -1;
  }
};

#endif // MATERIAL_SIGN_CAPSULE_H
