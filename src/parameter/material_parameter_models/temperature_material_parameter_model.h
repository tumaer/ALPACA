//===-------------- temperature_material_parameter_model.h ----------------===//
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
#ifndef TEMPERATURE_MATERIAL_PARAMETER_MODEL_H
#define TEMPERATURE_MATERIAL_PARAMETER_MODEL_H

#include "parameter/material_parameter_model.h"

/**
 * @brief The TemperatureMaterialParameterModel class implements the specific
 * computation for a material parameter model that acts on the prime state
 * temperature. It provides the loop structure to compute the parameter for a
 * single given block. The actual implementation of the model param = f(T) is
 * done in the derived class.
 *
 * @tparam DerivedTemperatureMaterialParameterModel derived temperature model
 * class that provides the model computation.
 */
template <typename DerivedTemperatureMaterialParameterModel>
class TemperatureMaterialParameterModel : public MaterialParameterModel {

  // friend model, which effectively implements the computation of the parameter
  friend DerivedTemperatureMaterialParameterModel;

  /**
   * @brief Executes the actual parameter calculation on the complete block.
   * @param block Block on which the parameter calculation should be carried out
   * (parameter on block as indirect return).
   */
  void DoUpdateParameter(Block &block, double const) const override {

    // extract the temperature from the block
    double const(&T)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::Temperature);

    // extract the paramerer from the block, which should be computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetParameterBuffer(
            DerivedTemperatureMaterialParameterModel::parameter_to_calculate_);

    // Compute the parameter based on the temperature
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
          parameter_buffer[i][j][k] =
              static_cast<DerivedTemperatureMaterialParameterModel const &>(
                  *this)
                  .ComputeParameter(T[i][j][k]);
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

    // extract the temperature from the block
    double const(&T)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::Temperature);

    // extract the parameter from the block, which should be computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetParameterBuffer(
            DerivedTemperatureMaterialParameterModel::parameter_to_calculate_);

    // Compute the parameter based on the shear rate
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {

          // Note: Volume fraction can be used to average temperature
          parameter_buffer[i][j][k] =
              interface_tags[i][j][k] * material_sign >= 0
                  ? static_cast<
                        DerivedTemperatureMaterialParameterModel const &>(*this)
                        .ComputeParameter(T[i][j][k])
                  : 0.0;

        } // k
      }   // j
    }     // i
  }

  /**
   * @brief Protected constructor to create the material temperature model.
   * @note Can only be called from derived classes.
   */
  explicit TemperatureMaterialParameterModel() : MaterialParameterModel() {
    /** Empty besides base class constructor call */
  }

public:
  virtual ~TemperatureMaterialParameterModel() = default;
  TemperatureMaterialParameterModel(TemperatureMaterialParameterModel const &) =
      delete;
  TemperatureMaterialParameterModel &
  operator=(TemperatureMaterialParameterModel const &) = delete;
  TemperatureMaterialParameterModel(TemperatureMaterialParameterModel &&) =
      delete;
  TemperatureMaterialParameterModel &
  operator=(TemperatureMaterialParameterModel &&) = delete;
};

#endif // TEMPERATURE_MATERIAL_PARAMETER_MODEL_H
