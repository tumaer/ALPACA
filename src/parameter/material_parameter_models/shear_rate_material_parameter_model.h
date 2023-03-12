//===-------------- shear_rate_material_parameter_model.h -----------------===//
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
#ifndef SHEAR_RATE_MATERIAL_PARAMETER_MODEL_H
#define SHEAR_RATE_MATERIAL_PARAMETER_MODEL_H

#include "parameter/material_parameter_model.h"
#include "stencils/stencil_utilities.h"
#include "utilities/shear_rate_operations.h"

/**
 * @brief The ShearRateMaterialParameterModel class implements the specific
 * computation for a material parameter model that acts on the shear rate. It
 * provides the loop structure to compute the parameter for a single given
 * block. The actual implementation of the model param = f(gamma_dot) is done in
 * the derived class.
 *
 * @tparam DerivedShearRateMaterialParameterModel derived shear rate model class
 * that provides the model computation.
 */
template <typename DerivedShearRateMaterialParameterModel>
class ShearRateMaterialParameterModel : public MaterialParameterModel {

  // friend model, which effectively implements the computation of the parameter
  friend DerivedShearRateMaterialParameterModel;

  /**
   * @brief Executes the actual parameter calculation on the complete block.
   * @param block Block on which the parameter calculation should be carried out
   * (parameter on block as indirect return).
   * @param cell_size Cell_size of the given block.
   */
  void DoUpdateParameter(Block &block, double const cell_size) const override {

    // Define the correct type of the stencil
    using DerivativeStencil = typename DerivativeStencilSetup::Concretize<
        DerivedShearRateMaterialParameterModel::derivative_stencil_>::type;

    // extract the velocities from the block
    // y and z velocity buffers may not be available
    // as workaround use the x velocity buffer in these cases
    // this is legal since the respective gradients are not computed/used anyway
    double const(&u)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::VelocityX);
    double const(&v)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::VelocityY);
    double const(&w)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::VelocityZ);

    // extract the shear viscosity from the block, which should be computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetParameterBuffer(
            DerivedShearRateMaterialParameterModel::parameter_buffer_type_);

    /**
     * Description for the positions of the Array:
     * [CC::TCX()]    [CC::TCY()]    [CC::TCZ()]
     * Field index x  Field index y  Field index z
     */
    double shear_rate[CC::TCX()][CC::TCY()][CC::TCZ()];

    // Initialize the shear_rate buffer
    for (unsigned int i = 0; i < CC::TCX(); ++i) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int k = 0; k < CC::TCZ(); ++k) {

          shear_rate[i][j][k] = 0.0;
        }
      }
    }

    // Compute the shear rate from the velocities
    ShearRateOperations::ComputeShearRate<DerivativeStencil>(u, v, w, cell_size,
                                                             shear_rate);

    // Compute the parameter based on the shear rate
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
          parameter_buffer[i][j][k] =
              static_cast<DerivedShearRateMaterialParameterModel const &>(*this)
                  .ComputeParameter(shear_rate[i][j][k]);
        }
      }
    }
  }

  /**
   * @brief Executes the actual parameter calculation on the on the block for
   * the given material up to the interface.
   * @param block Block on which the parameter calculation should be carried out
   * (parameter on block as indirect return).
   * @param cell_size Cell_size of the given block.
   * @param interface_tags Tags describing the interface position.
   * @param material_sign Sign of the material for identification on interface
   * tags.
   */
  void DoUpdateParameter(
      Block &block, double const cell_size,
      std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
      std::int8_t const material_sign) const override {

    // Define the correct type of the stencil
    using DerivativeStencil = typename DerivativeStencilSetup::Concretize<
        DerivedShearRateMaterialParameterModel::derivative_stencil_>::type;

    // extract the velocities from the block
    double const(&u)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::VelocityX);
    double const(&v)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::VelocityY);
    double const(&w)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetPrimeStateBuffer(PrimeState::VelocityZ);

    // extract the appropriate parameter buffer from the block, which should be
    // computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        block.GetParameterBuffer(
            DerivedShearRateMaterialParameterModel::parameter_buffer_type_);

    // Compute the viscosity based on the shear rate and model parameter
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
          if (interface_tags[i][j][k] * material_sign >= 0) {
            // compute the shear rate for this cell
            double const shear_rate =
                ShearRateOperations::ComputeShearRate<DerivativeStencil>(
                    u, v, w, cell_size, i, j, k);
            // compute the parameter for the cell
            parameter_buffer[i][j][k] =
                static_cast<DerivedShearRateMaterialParameterModel const &>(
                    *this)
                    .ComputeParameter(shear_rate);
          } else {
            parameter_buffer[i][j][k] = 0.0;
          }
        } // k
      }   // j
    }     // i
  }

  /**
   * @brief Protected constructor to create the material shear rate model.
   * @note Can only be called from derived classes.
   */
  explicit ShearRateMaterialParameterModel() : MaterialParameterModel() {
    /** Empty besides base class constructor call */
  }

public:
  virtual ~ShearRateMaterialParameterModel() = default;
  ShearRateMaterialParameterModel(ShearRateMaterialParameterModel const &) =
      delete;
  ShearRateMaterialParameterModel &
  operator=(ShearRateMaterialParameterModel const &) = delete;
  ShearRateMaterialParameterModel(ShearRateMaterialParameterModel &&) = delete;
  ShearRateMaterialParameterModel &
  operator=(ShearRateMaterialParameterModel &&) = delete;
};

#endif // SHEAR_RATE_MATERIAL_PARAMETER_MODEL_H
