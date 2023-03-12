//===----------- constant_surface_tension_coefficient_model.h -------------===//
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
#ifndef CONSTANT_SURFACE_TENSION_COEFFICIENT_MODEL_H
#define CONSTANT_SURFACE_TENSION_COEFFICIENT_MODEL_H

#include "parameter/interface_parameter_models/constant_interface_parameter_model.h"
#include "unit_handler.h"

/**
 * @brief Class that provides the actual model implementation f = const. for a
 * Constant material parameter model to compute the thermal conductivity in the
 * field.
 */
class ConstantSurfaceTensionCoefficientModel
    : public ConstantInterfaceParameterModel<
          ConstantSurfaceTensionCoefficientModel> {

  // Defines the friend to be used for the computation indicating the
  // primestates this model depends on
  friend ConstantInterfaceParameterModel;

  // Definition of parameters needed in the base class for the computations
  static constexpr InterfaceParameter parameter_buffer_type_ =
      InterfaceParameter::SurfaceTensionCoefficient;

  // member variables needed for the model calculation
  double const sigma_constant_;

  // Actual computation of surface tension coefficient
  double ComputeParameter() const;

public:
  ConstantSurfaceTensionCoefficientModel() = delete;
  explicit ConstantSurfaceTensionCoefficientModel(
      std::unordered_map<std::string, double> const &parameter_map,
      UnitHandler const &unit_handler);
  virtual ~ConstantSurfaceTensionCoefficientModel() = default;
  ConstantSurfaceTensionCoefficientModel(
      const ConstantSurfaceTensionCoefficientModel &) = delete;
  ConstantSurfaceTensionCoefficientModel &
  operator=(const ConstantSurfaceTensionCoefficientModel &) = delete;
  ConstantSurfaceTensionCoefficientModel(
      ConstantSurfaceTensionCoefficientModel &&) = delete;
  ConstantSurfaceTensionCoefficientModel &
  operator=(ConstantSurfaceTensionCoefficientModel &&) = delete;

  // logging function
  std::string GetLogData(unsigned int const indent,
                         UnitHandler const &unit_handler) const;
};

#endif // CONSTANT_SURFACE_TENSION_COEFFICIENT_MODEL_H
