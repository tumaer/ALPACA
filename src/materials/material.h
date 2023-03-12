//===--------------------------- material.h -------------------------------===//
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
#ifndef MATERIAL_H
#define MATERIAL_H

#include <memory>
#include <vector>

#include "materials/equation_of_state.h"
#include "materials/material_definitions.h"
#include "parameter/material_parameter_model.h"
#include "unit_handler.h"

/**
 * @brief The Material class defines an interface for different materials
 * containing an equation of state and other pure material properties (e.g.
 * viscosity, thermal conductivity, specific heat capacity). Furthermore, it
 * stores parameter models for properties that allow such computations. The
 * Material class does not manipulate any data.
 */
class Material {

  // const pointer to the specific equation of state
  std::unique_ptr<EquationOfState const> equation_of_state_;

  // material properties (fixed values).
  double const bulk_viscosity_;
  double const shear_viscosity_;
  double const thermal_conductivity_;
  double const specific_heat_capacity_;
  // material property models (non-fixed computations)
  std::unique_ptr<MaterialParameterModel const> shear_viscosity_model_;
  std::unique_ptr<MaterialParameterModel const> thermal_conductivity_model_;

public:
  explicit Material(
      std::unique_ptr<EquationOfState const> equation_of_state,
      double const dimensional_bulk_viscosity,
      double const dimensional_shear_viscosity,
      double const dimensional_thermal_conductivity,
      double const dimensional_specific_heat_capacity,
      std::unique_ptr<MaterialParameterModel const> shear_viscosity_model,
      std::unique_ptr<MaterialParameterModel const> thermal_conductivity_model,
      UnitHandler const &unit_handler);
  Material() = delete;
  ~Material() = default;
  Material(Material const &) = delete;
  Material &operator=(Material const &) = delete;
  // Non-deleted move constructor
  Material(Material &&material);
  Material &operator=(Material &&) = delete;

  // return functions of member variables
  EquationOfState const &GetEquationOfState() const;
  std::vector<double> GetShearAndBulkViscosity() const;
  double GetShearViscosity() const;
  double GetBulkViscosity() const;
  double GetThermalConductivity() const;
  double GetSpecificHeatCapacity() const;
  MaterialParameterModel const &GetShearViscosityModel() const;
  MaterialParameterModel const &GetThermalConductivityModel() const;
};

#endif
