//===----------------- power_law_shear_viscosity_model.h ------------------===//
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
#ifndef POWER_LAW_SHEAR_VISCOSITY_MODEL_H
#define POWER_LAW_SHEAR_VISCOSITY_MODEL_H

#include "parameter/material_parameter_models/shear_rate_material_parameter_model.h"
#include "unit_handler.h"

/**
 * @brief Class that provides the actual model implementation param =
 * f(gamma_dot) following a power law. for a shear-rate based material parameter
 * model to compute the shear viscosity in the field. The model is implemented
 * based on \cite{Deville2012}.
 */
class PowerLawShearViscosityModel
    : public ShearRateMaterialParameterModel<PowerLawShearViscosityModel> {

  // Defines the friend to be used for the computation indicating the variables
  // this model depends on (e.g., prime states)
  friend ShearRateMaterialParameterModel;

  // Definition of parameters needed in the base class for the computations
  static constexpr Parameter parameter_buffer_type_ =
      Parameter::ShearViscosity; // buffer that is modified by the model
  static constexpr DerivativeStencils derivative_stencil_ =
      viscous_fluxes_derivative_stencil_cell_center; // stencil required for the
                                                     // shear rate

  // member variables needed for the model calculation
  double const consistency_factor_;
  double const power_law_exponent_;
  // derived model parameter (to avoid multiple computations)
  double const exponent_;

  // Actual computation of shear viscosity
  double ComputeParameter(double const shear_rate) const;

public:
  PowerLawShearViscosityModel() = delete;
  explicit PowerLawShearViscosityModel(
      std::unordered_map<std::string, double> const &dimensional_parameter_map,
      UnitHandler const &unit_handler);
  virtual ~PowerLawShearViscosityModel() = default;
  PowerLawShearViscosityModel(const PowerLawShearViscosityModel &) = delete;
  PowerLawShearViscosityModel &
  operator=(const PowerLawShearViscosityModel &) = delete;
  PowerLawShearViscosityModel(PowerLawShearViscosityModel &&) = delete;
  PowerLawShearViscosityModel &
  operator=(PowerLawShearViscosityModel &&) = delete;

  // logging function
  std::string GetLogData(unsigned int const indent,
                         UnitHandler const &unit_handler) const;
};

#endif // POWER_LAW_SHEAR_VISCOSITY_MODEL_H
