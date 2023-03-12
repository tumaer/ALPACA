//===---------------------- material_pairing.h ----------------------------===//
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
#ifndef MATERIAL_PAIRING_H
#define MATERIAL_PAIRING_H

#include "parameter/interface_parameter_model.h"
#include "unit_handler.h"

/**
 * @brief The MaterialPairing class defines an interface for all pairing
 * properties used in the simulation (e.g. surface tension coefficient).
 *        Furthermore, it stores parameter models for properties that allow such
 * computations. The MaterialPairing class does not manipulate any data.
 */
class MaterialPairing {

protected:
  // material pairing properties (fixed values).
  double const surface_tension_coefficient_;
  // models
  std::unique_ptr<InterfaceParameterModel const>
      surface_tension_coefficient_model_;

public:
  explicit MaterialPairing(double const dimensional_surface_tension_coefficient,
                           std::unique_ptr<InterfaceParameterModel const>
                               surface_tension_coefficient_model,
                           UnitHandler const &unit_handler);
  explicit MaterialPairing();
  virtual ~MaterialPairing() = default;
  MaterialPairing(MaterialPairing const &) = delete;
  MaterialPairing &operator=(MaterialPairing const &) = delete;
  // Non-deleted move constructor
  MaterialPairing(MaterialPairing &&);
  MaterialPairing &operator=(MaterialPairing &&) = delete;

  // return functions of member variables
  double GetSurfaceTensionCoefficient() const;
  InterfaceParameterModel const &GetSurfaceTensionCoefficientModel() const;
};

#endif // MATERIAL_PAIRING_H
