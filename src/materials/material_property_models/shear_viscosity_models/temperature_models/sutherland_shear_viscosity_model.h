//===--------------- sutherland_shear_viscosity_model.h -------------------===//
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
#ifndef SUTHERLAND_SHEAR_SHEAR_VISCOSITY_MODEL_H
#define SUTHERLAND_SHEAR_SHEAR_VISCOSITY_MODEL_H

#include "block_definitions/field_material_definitions.h"
#include "parameter/material_parameter_models/temperature_material_parameter_model.h"
#include "unit_handler.h"
#include <cmath>

/**
 * @brief Class that provides the actual model implementation param = f(T)
 * following Sutherlands viscosity law. for a temperature based material
 * parameter model to compute the shear viscosity in the field. The model is
 * implemented base on \cite{Sutherland1893}.
 */
class SutherlandShearViscosityModel
    : public TemperatureMaterialParameterModel<SutherlandShearViscosityModel> {

  // Defines the friend to be used for the computation indicating the variables
  // this model depends on (e.g., prime states)
  friend TemperatureMaterialParameterModel;

  // Definition of parameters needed in the base class for the computations
  static constexpr Parameter parameter_to_calculate_ =
      Parameter::ShearViscosity;

  // member variables needed for the model calculation
  double const mu0_;
  double const T0_;
  double const S_;
  // derived model parameter (to avoid multiple computations)
  double const constant_;
  double const epsilon_ = std::numeric_limits<double>::epsilon();

  // Actual computation of shear viscosity
  double ComputeParameter(double const temperature) const;

public:
  SutherlandShearViscosityModel() = delete;
  explicit SutherlandShearViscosityModel(
      std::unordered_map<std::string, double> const &dimensional_parameter_map,
      UnitHandler const &unit_handler);
  virtual ~SutherlandShearViscosityModel() = default;
  SutherlandShearViscosityModel(const SutherlandShearViscosityModel &) = delete;
  SutherlandShearViscosityModel &
  operator=(const SutherlandShearViscosityModel &) = delete;
  SutherlandShearViscosityModel(SutherlandShearViscosityModel &&) = delete;
  SutherlandShearViscosityModel &
  operator=(SutherlandShearViscosityModel &&) = delete;

  // logging function
  std::string GetLogData(unsigned int const indent,
                         UnitHandler const &unit_handler) const;
};

#endif // SUTHERLAND_SHEAR_SHEAR_VISCOSITY_MODEL_H
