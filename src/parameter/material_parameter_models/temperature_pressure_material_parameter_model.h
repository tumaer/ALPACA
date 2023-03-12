//===---------- temperature_pressure_material_parameter_model.h -----------===//
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
#ifndef TEMPERATURE_PRESSURE_MATERIAL_PARAMETER_MODEL_H
#define TEMPERATURE_PRESSURE_MATERIAL_PARAMETER_MODEL_H

#include "parameter/material_parameter_model.h"

/**
 * @brief The TemperaturePressureMaterialParameterModel class implements the
 * specific computation for a material parameter model that acts on the prime
 * states temperature and pressure. It provides the loop structure to compute
 * the parameter for a single given block. The actual implementation of the
 *        model param = f(T,p) is done in the derived class.
 *
 * @tparam DerivedTemperaturePressureMaterialParameterModel derived temperature
 * pressure model class that provides the model computations.
 */
template <typename DerivedTemperaturePressureMaterialParameterModel>
class TemperaturePressureMaterialParameterModel
    : public MaterialParameterModel {

  // friend model, which effectively implements the computation of the parameter
  friend DerivedTemperaturePressureMaterialParameterModel;

  /**
   * @brief Executes the actual parameter calculation on the complete block.
   * @param block Block on which the parameter calculation should be carried out
   * (parameter on block as indirect return).
   */
  void DoUpdateParameter(Block &block, double const) const override {

    // extract the temperature from the block
    double const(&T)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::Temperature);
    double const(&p)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::Pressure);

    // extract the shear viscosity from the block, which should be computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetParameterBuffer(
            DerivedTemperaturePressureMaterialParameterModel::
                parameter_to_calculate_);

    // Compute the parameter based on the temperature
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {

          parameter_buffer[i][j][k] =
              static_cast<
                  DerivedTemperaturePressureMaterialParameterModel const &>(
                  *this)
                  .ComputeParameter(T[i][j][k], p[i][j][k]);
        }
      }
    }
  }

  /**
   * @brief Executes the actual parameter calculation on the on the block for
   * the given material up to the interface.
   * @param block Block on which the parameter calculation should be carried out
   * (parameter on block as indirect return).
   * @param interface_tags Tags describing the interface position.
   * @param material_sign Sign of the material for identification on interface
   * tags.
   */
  void DoUpdateParameter(
      Block &block, double const,
      std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
      std::int8_t const material_sign) const override {

    // extract the temperature and pressure from the block
    double const(&T)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::Temperature);
    double const(&p)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::Pressure);

    // extract the parameter from the block, which should be computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetParameterBuffer(
            DerivedTemperaturePressureMaterialParameterModel::
                parameter_to_calculate_);

    // Compute the parameter based on the shear rate and model parameter
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {

          // Note: Volume fraction can be used to average temperature and
          // pressure
          parameter_buffer[i][j][k] =
              interface_tags[i][j][k] * material_sign >= 0
                  ? static_cast<
                        DerivedTemperaturePressureMaterialParameterModel const
                            &>(*this)
                        .ComputeParameter(T[i][j][k], p[i][j][k])
                  : 0.0;

        } // k
      }   // j
    }     // i
  }

  /**
   * @brief Protected constructor to create the material temperature pressure
   * model.
   * @note Can only be called from derived classes.
   */
  explicit TemperaturePressureMaterialParameterModel()
      : MaterialParameterModel() {
    /** Empty besides base class constructor call */
  }

public:
  virtual ~TemperaturePressureMaterialParameterModel() = default;
  TemperaturePressureMaterialParameterModel(
      TemperaturePressureMaterialParameterModel const &) = delete;
  TemperaturePressureMaterialParameterModel &
  operator=(TemperaturePressureMaterialParameterModel const &) = delete;
  TemperaturePressureMaterialParameterModel(
      TemperaturePressureMaterialParameterModel &&) = delete;
  TemperaturePressureMaterialParameterModel &
  operator=(TemperaturePressureMaterialParameterModel &&) = delete;
};

#endif // TEMPERATURE_PRESSURE_MATERIAL_PARAMETER_MODEL_H
