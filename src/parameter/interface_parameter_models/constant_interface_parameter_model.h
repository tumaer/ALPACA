//===--------------- constant_interface_parameter_model.h -----------------===//
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
#ifndef CONSTANT_INTERFACE_PARAMETER_MODEL_H
#define CONSTANT_INTERFACE_PARAMETER_MODEL_H

#include "enums/interface_tag_definition.h"
#include "parameter/interface_parameter_model.h"

/**
 * @brief The ConstantInterfaceParameterModel class implements the specific
 * computation for a interface parameter model that gives constant values. It
 * provides the loop structure to compute the parameter for a single given node.
 * The actual implementation of the model param = const. is done in the derived
 * class.
 * @note This class serves as a compatability class to provide the same behavior
 * for all material pairings in case one uses a parameter model to avoid
 * additional if else statements calling if a model is used or not. Therefore,
 * if one material uses a model, all materials at least must define a constant
 * model.
 *
 * @tparam DerivedConstantInterfaceParameterModel derived constan model class
 * that provides the model computation.
 */
template <typename DerivedConstantInterfaceParameterModel>
class ConstantInterfaceParameterModel : public InterfaceParameterModel {

  // friend model, which effectively implements the computation of the parameter
  friend DerivedConstantInterfaceParameterModel;

  /**
   * @brief Executes the actual parameter calculation on the complete block
   * @param node Node for which the parameter is computed
   */
  void DoUpdateParameter(Node &node) const override {

    // Obtain the interface tags and volume fraction
    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

    // extract the parameter from the block, which should be computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetInterfaceParameterBuffer(
            DerivedConstantInterfaceParameterModel::parameter_buffer_type_);

    // Compute the parameter based on constan values
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {

          // Only compute the parameter for cut cells. Average quantities based
          // on the volume fraction
          if (interface_tags[i][j][k] == ITTI(IT::OldCutCell)) {
            parameter_buffer[i][j][k] =
                static_cast<DerivedConstantInterfaceParameterModel const &>(
                    *this)
                    .ComputeParameter();
          } else {
            parameter_buffer[i][j][k] = 0.0;
          }
        } // k
      }   // j
    }     // i
  }

  /**
   * @brief Protected constructor to create the interface constant model.
   * @note Can only be called from derived classes.
   */
  explicit ConstantInterfaceParameterModel() : InterfaceParameterModel() {
    /** Empty besides base class constructor call */
  }

public:
  virtual ~ConstantInterfaceParameterModel() = default;
  ConstantInterfaceParameterModel(ConstantInterfaceParameterModel const &) =
      delete;
  ConstantInterfaceParameterModel &
  operator=(ConstantInterfaceParameterModel const &) = delete;
  ConstantInterfaceParameterModel(ConstantInterfaceParameterModel &&) = delete;
  ConstantInterfaceParameterModel &
  operator=(ConstantInterfaceParameterModel &&) = delete;
};

#endif // CONSTANT_INTERFACE_PARAMETER_MODEL_H
