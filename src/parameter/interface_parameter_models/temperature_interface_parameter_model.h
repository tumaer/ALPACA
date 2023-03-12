//===----------- temperature_interface_parameter_model.h ------------------===//
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
#ifndef TEMPERATURE_INTERFACE_PARAMETER_MODEL_H
#define TEMPERATURE_INTERFACE_PARAMETER_MODEL_H

#include "parameter/interface_parameter_model.h"

/**
 * @brief The TemperatureInterfaceParameterModel class implements the specific
 * computation for a interface parameter model that gives temperature dependent
 * values. It provides the loop structure to compute the parameter for a single
 * given node. The actual implementation of the model param = f(T) is done in
 * the derived class.
 *
 * @tparam DerivedTemperatureInterfaceParameterModel derived temperature model
 * class that provides the model computation.
 */
template <typename DerivedTemperatureInterfaceParameterModel>
class TemperatureInterfaceParameterModel : public InterfaceParameterModel {

  // friend model, which effectively implements the computation of the parameter
  friend DerivedTemperatureInterfaceParameterModel;

  /**
   * @brief Executes the actual parameter calculation on the complete block.
   * @param node Node for which the parameter is computed.
   */
  void DoUpdateParameter(Node &node) const override {

    // extract the temperature from the block
    double const(&T_positive)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetPhaseByMaterial(MaterialSignCapsule::PositiveMaterial())
            .GetPrimeStateBuffer(PrimeState::Temperature);
    double const(&T_negative)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetPhaseByMaterial(MaterialSignCapsule::NegativeMaterial())
            .GetPrimeStateBuffer(PrimeState::Temperature);

    // extract the parameter from the block, which should be computed
    double(&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetInterfaceParameterBuffer(
            DerivedTemperatureInterfaceParameterModel::parameter_buffer_type_);

    // Obtain the interface tags and volume fraction
    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
    double const(&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetBaseBuffer(Levelset::VolumeFraction);

    // Compute the parameter based on constan values
    for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {

          // Only compute the parameter for cut cells. Average quantities based
          // on the volume fraction
          if (interface_tags[i][j][k] == ITTI(IT::OldCutCell)) {
            double const averaged_T =
                volume_fraction[i][j][k] * T_positive[i][j][k] +
                (1.0 - volume_fraction[i][j][k]) * T_negative[i][j][k];
            parameter_buffer[i][j][k] =
                static_cast<DerivedTemperatureInterfaceParameterModel const &>(
                    *this)
                    .ComputeParameter(averaged_T);
          } else {
            parameter_buffer[i][j][k] = 0.0;
          }
        } // k
      }   // j
    }     // i
  }

  /**
   * @brief Protected constructor to create the interface temperature model.
   * @note Can only be called from derived classes.
   */
  explicit TemperatureInterfaceParameterModel() : InterfaceParameterModel() {
    /** Empty besides base class constructor call */
  }

public:
  virtual ~TemperatureInterfaceParameterModel() = default;
  TemperatureInterfaceParameterModel(
      TemperatureInterfaceParameterModel const &) = delete;
  TemperatureInterfaceParameterModel &
  operator=(TemperatureInterfaceParameterModel const &) = delete;
  TemperatureInterfaceParameterModel(TemperatureInterfaceParameterModel &&) =
      delete;
  TemperatureInterfaceParameterModel &
  operator=(TemperatureInterfaceParameterModel &&) = delete;
};

#endif // TEMPERATURE_INTERFACE_PARAMETER_MODEL_H
