//===---------------- constant_shear_viscosity_model.h --------------------===//
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
#ifndef CONSTANT_SHEAR_VISCOSITY_MODEL_H
#define CONSTANT_SHEAR_VISCOSITY_MODEL_H

#include "parameter/material_parameter_models/constant_material_parameter_model.h"
#include "unit_handler.h"

/**
 * @brief Class that provides the actual model implementation f = const. for a
 * constant material parameter model to compute the field (dependent) shear
 * viscosity. This model only serves as a complement to the other material
 * parameter models to provide field dependent computations for all materials if
 * at least one uses a real parameter model.
 */
class ConstantShearViscosityModel
    : public ConstantMaterialParameterModel<ConstantShearViscosityModel> {

  // Defines the friend to be used for the computation indicating the variables
  // this model depends on (e.g., prime states)
  friend ConstantMaterialParameterModel;

  // Definition of all parameters needed in the base class for the computations
  static constexpr Parameter parameter_buffer_type_ =
      Parameter::ShearViscosity; // buffer that is modified by the model

  // member variables needed for the model calculation
  double const mu_constant_;

  // Actual computation of shear viscosity
  double ComputeParameter() const;

public:
  ConstantShearViscosityModel() = delete;
  explicit ConstantShearViscosityModel(
      std::unordered_map<std::string, double> const &dimensional_parameter_map,
      UnitHandler const &unit_handler);
  virtual ~ConstantShearViscosityModel() = default;
  ConstantShearViscosityModel(const ConstantShearViscosityModel &) = delete;
  ConstantShearViscosityModel &
  operator=(const ConstantShearViscosityModel &) = delete;
  ConstantShearViscosityModel(ConstantShearViscosityModel &&) = delete;
  ConstantShearViscosityModel &
  operator=(ConstantShearViscosityModel &&) = delete;

  // logging function
  std::string GetLogData(unsigned int const indent,
                         UnitHandler const &unit_handler) const;
};

#endif // CONSTANT_SHEAR_VISCOSITY_MODEL_H
