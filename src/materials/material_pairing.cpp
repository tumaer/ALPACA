//===---------------------- material_pairing.cpp --------------------------===//
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
#include "materials/material_pairing.h"

#include "materials/material_property_definitions.h"
#include "utilities/string_operations.h"

/**
 * @brief Sets up a material pairing for the given input data.
 * @param dimensional_surface_tension_coefficient Fixed dimensional value of the
 * surface tension coefficient (ownership takes place).
 * @param surface_tension_coefficient_model The model used for the surface
 * tension coefficient.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 *
 * @note No default values are specified to ensure that everything is provided.
 * Default values are set during the initialization.
 */
MaterialPairing::MaterialPairing(
    double const dimensional_surface_tension_coefficient,
    std::unique_ptr<InterfaceParameterModel const>
        surface_tension_coefficient_model,
    UnitHandler const &unit_handler)
    : // Start initializer list
      surface_tension_coefficient_(unit_handler.NonDimensionalizeValue(
          dimensional_surface_tension_coefficient,
          UnitType::SurfaceTensionCoefficient)),
      surface_tension_coefficient_model_(
          std::move(surface_tension_coefficient_model)) {
  /** Empty besides initializer list */
}

/**
 * @brief Sets up a material pairing with no input data
 * @note Function can be removed when the surface tension coefficient hard
 * coding is removed from the capillary pressure calculator.
 */
MaterialPairing::MaterialPairing()
    : surface_tension_coefficient_(0.0),
      surface_tension_coefficient_model_(nullptr) {
  /** Empty besides initializer list */
}

/**
 * @brief Move constructor.
 * @param pairing Material pairing from which the data are moved.
 */
MaterialPairing::MaterialPairing(MaterialPairing &&pairing)
    : // Start initializer list
      surface_tension_coefficient_(pairing.surface_tension_coefficient_),
      surface_tension_coefficient_model_(
          std::move(pairing.surface_tension_coefficient_model_)) {
  /** Empty besides initializer list */
}

/**
 * @brief Gives the surface tension coefficient  (fixed value)
 * @return surface tension coefficient
 */
double MaterialPairing::GetSurfaceTensionCoefficient() const {
  return surface_tension_coefficient_;
}

/**
 * @brief Gives the model of the surface tension coefficient.
 * @return The instance of the surface tension coefficient model.
 */
InterfaceParameterModel const &
MaterialPairing::GetSurfaceTensionCoefficientModel() const {
  return *surface_tension_coefficient_model_;
}
