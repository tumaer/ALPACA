//===--------------- constant_material_parameter_model.h ------------------===//
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
#ifndef CONSTANT_MATERIAL_PARAMETER_MODEL_H
#define CONSTANT_MATERIAL_PARAMETER_MODEL_H

#include "parameter/material_parameter_model.h"

/**
 * @brief The ConstantMaterialParameterModel class implements the specific
 * computation for a material parameter model that gives constant values. It
 * provides the loop structure to compute the parameter for a single given
 * block. The actual implementation of the model param = const. is done in the
 * derived class.
 * @note This class serves as a compatability class to provide the same behavior
 * for all materials in case one uses a parameter model to avoid additional if
 * else statements calling if a model is used or not. Therefore, if one material
 * uses a model, all materials at least must define a constant model.
 *
 * @tparam DerivedConstantMaterialParameterModel derived constan model class
 * that provides the model computation.
 */
template <typename DerivedConstantMaterialParameterModel>
class ConstantMaterialParameterModel : public MaterialParameterModel {

  // friend model, which effectively implements the computation of the parameter
  friend DerivedConstantMaterialParameterModel;

  /**
   * @brief Executes the actual parameter calculation on the complete block.
   * @param block Block on which the parameter calculation should be carried out
   * (parameter on block as indirect return).
   */
  void DoUpdateParameter(Block &block, double const) const override {

    // extract the parameter from the block, which should be computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetParameterBuffer(
            DerivedConstantMaterialParameterModel::parameter_buffer_type_);

    // Compute the parameter based on constant values
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
          parameter_buffer[i][j][k] =
              static_cast<DerivedConstantMaterialParameterModel const &>(*this)
                  .ComputeParameter();
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

    // extract the parameter from the block, which should be computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetParameterBuffer(
            DerivedConstantMaterialParameterModel::parameter_buffer_type_);

    // Compute the parameter based on constan values
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {

          parameter_buffer[i][j][k] =
              interface_tags[i][j][k] * material_sign >= 0
                  ? static_cast<DerivedConstantMaterialParameterModel const &>(
                        *this)
                        .ComputeParameter()
                  : 0.0;
        } // k
      }   // j
    }     // i
  }

  /**
   * @brief Protected constructor to create the material constant model.
   * @note Can only be called from derived classes.
   */
  explicit ConstantMaterialParameterModel() : MaterialParameterModel() {
    /** Empty besides base class constructor call */
  }

public:
  virtual ~ConstantMaterialParameterModel() = default;
  ConstantMaterialParameterModel(ConstantMaterialParameterModel const &) =
      delete;
  ConstantMaterialParameterModel &
  operator=(ConstantMaterialParameterModel const &) = delete;
  ConstantMaterialParameterModel(ConstantMaterialParameterModel &&) = delete;
  ConstantMaterialParameterModel &
  operator=(ConstantMaterialParameterModel &&) = delete;
};

#endif // CONSTANT_MATERIAL_PARAMETER_MODEL_H
