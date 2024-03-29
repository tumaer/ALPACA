//===--------------------- prime_state_handler.h --------------------------===//
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
#ifndef PRIME_STATE_HANDLER_H
#define PRIME_STATE_HANDLER_H

#include "block_definitions/block.h"
#include "block_definitions/field_material_definitions.h"
#include "materials/equations_of_state/gamma_model_stiffened_gas.h"
#include "materials/material_manager.h"
#include "user_specifications/equation_settings.h"

/**
 * @brief This class handles the conversion between conservative and prime state
 * values for full buffers or single cells.
 */
class PrimeStateHandler {

  MaterialManager const &material_manager_;

  template <typename PrimeStatesContainerType,
            typename ConservativesContainerType>
  inline void DoConvertPrimeStatesToConservatives(
      MaterialName const &material,
      PrimeStatesContainerType const &prime_states_container,
      ConservativesContainerType &conservatives_container) const {

    PrimeStatesToConservatives(
        material_manager_.GetMaterial(material).GetEquationOfState(),
        prime_states_container, conservatives_container);
  }

  template <typename ConservativesContainerType,
            typename PrimeStatesContainerType>
  inline void DoConvertConservativesToPrimeStates(
      MaterialName const &material,
      ConservativesContainerType const &conservatives_container,
      PrimeStatesContainerType &prime_states_container) const {

    ConservativesToPrimeStates(
        material_manager_.GetMaterial(material).GetEquationOfState(),
        conservatives_container, prime_states_container);
  }

public:
  explicit PrimeStateHandler(MaterialManager const &material_manager)
      : material_manager_(material_manager) {}
  ~PrimeStateHandler() = default;
  PrimeStateHandler(PrimeStateHandler const &) = delete;
  PrimeStateHandler &operator=(PrimeStateHandler const &) = delete;
  PrimeStateHandler(PrimeStateHandler &&) = delete;
  PrimeStateHandler &operator=(PrimeStateHandler &&) = delete;

  /**
   * @brief Converts prime state to conservative values for every cell in the
   * given buffers.
   * @param material Material identifier of the material under consideration.
   * @param prime_states The input PrimeStates buffer.
   * @param conservatives The output Conservatives buffer.
   */
  void ConvertPrimeStatesToConservatives(MaterialName const &material,
                                         PrimeStates const &prime_states,
                                         Conservatives &conservatives) const {
    for (unsigned int i = 0; i < CC::TCX(); ++i) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int k = 0; k < CC::TCZ(); ++k) {
          ConvertPrimeStatesToConservatives(material, prime_states,
                                            conservatives, i, j, k);
        }
      }
    }
  }

  /**
   * @brief Converts prime state to conservative values for a specific cell in
   * the given buffers.
   * @param material Material identifier of the material under consideration.
   * @param prime_states The input PrimeStates buffer.
   * @param conservatives The output Conservatives buffer.
   * @param i,j,k Indices of the cell.
   */
  inline void ConvertPrimeStatesToConservatives(MaterialName const &material,
                                                PrimeStates const &prime_states,
                                                Conservatives &conservatives,
                                                unsigned int const i,
                                                unsigned int const j,
                                                unsigned int const k) const {
    auto &&conservatives_cell = conservatives.GetCellView(i, j, k);
    ConvertPrimeStatesToConservatives(
        material, prime_states.GetCellView(i, j, k), conservatives_cell);
  }

  /**
   * @brief Converts prime state to conservative values.
   * @param material Material identifier of the material under consideration.
   * @param prime_states_container Container of the input prime states.
   * @param conservatives_container Container of the output conservatives.
   * @tparam PrimeStatesContainerType Type of the conservative container, has to
   * provide an [index]-operator.
   * @tparam ConservativesContainerType Type of the prime state container, has
   * to provide an [index]-operator.
   * @note In the containers, prime state and conservative values are assumed to
   * be at the indices given by the PrimeState and Equation enumerations,
   * respectively.
   */
  template <typename PrimeStatesContainerType,
            typename ConservativesContainerType>
  inline void ConvertPrimeStatesToConservatives(
      MaterialName const &material,
      PrimeStatesContainerType const &prime_states_container,
      ConservativesContainerType &conservatives_container) const {
    DoConvertPrimeStatesToConservatives(material, prime_states_container,
                                        conservatives_container);
  }

  /**
   * @brief Converts conservative to prime state values for every cell in the
   * given buffers.
   * @param material Material identifier of the material under consideration.
   * @param conservatives The input Conservatives buffer.
   * @param prime_states The output PrimeStates buffer.
   */
  void ConvertConservativesToPrimeStates(MaterialName const &material,
                                         Conservatives const &conservatives,
                                         PrimeStates &prime_states) const {
    for (unsigned int i = 0; i < CC::TCX(); ++i) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int k = 0; k < CC::TCZ(); ++k) {
          ConvertConservativesToPrimeStates(material, conservatives,
                                            prime_states, i, j, k);
        }
      }
    }
  }

  /**
   * @brief Converts conservative to prime state values for a specific cell in
   * the given buffers.
   * @param material Material identifier of the material under consideration.
   * @param conservatives The input Conservatives buffer.
   * @param prime_states The output PrimeStates buffer.
   * @param i,j,k Indices of the cell.
   */
  inline void ConvertConservativesToPrimeStates(
      MaterialName const &material, Conservatives const &conservatives,
      PrimeStates &prime_states, unsigned int const i, unsigned int const j,
      unsigned int const k) const {
    auto &&prime_states_cell = prime_states.GetCellView(i, j, k);
    ConvertConservativesToPrimeStates(
        material, conservatives.GetCellView(i, j, k), prime_states_cell);
  }

  /**
   * @brief Converts conservative to prime state values.
   * @param material Material identifier of the material under consideration.
   * @param conservatives_container Container of the input conservatives.
   * @param prime_states_container Container of the output prime states.
   * @tparam ConservativesContainerType Type of the prime state container, has
   * to provide an [index]-operator.
   * @tparam PrimeStatesContainerType Type of the conservative container, has to
   * provide an [index]-operator.
   * @note In the containers, prime state and conservative values are assumed to
   * be at the indices given by the PrimeState and Equation enumerations,
   * respectively.
   */
  template <typename ConservativesContainerType,
            typename PrimeStatesContainerType>
  inline void ConvertConservativesToPrimeStates(
      MaterialName const &material,
      ConservativesContainerType const &conservatives_container,
      PrimeStatesContainerType &prime_states_container) const {
    DoConvertConservativesToPrimeStates(material, conservatives_container,
                                        prime_states_container);
  }
};

template <typename PrimeStatesContainerType,
          typename ConservativesContainerType>
inline void PrimeStatesToConservatives(
    EquationOfState const &eos,
    PrimeStatesContainerType const &prime_states_container,
    ConservativesContainerType &conservatives_container) {
  if constexpr (MF::IsEquationActive(Equation::Mass) &&
                MF::IsPrimeStateActive(PrimeState::Density)) {
    conservatives_container[ETI(Equation::Mass)] =
        prime_states_container[PTI(PrimeState::Density)];
  }
  if constexpr (MF::IsEquationActive(Equation::MomentumX) &&
                MF::IsPrimeStateActive(PrimeState::Density) &&
                MF::IsPrimeStateActive(PrimeState::VelocityX)) {
    conservatives_container[ETI(Equation::MomentumX)] =
        prime_states_container[PTI(PrimeState::VelocityX)] *
        prime_states_container[PTI(PrimeState::Density)];
  }
  if constexpr (MF::IsEquationActive(Equation::MomentumY) &&
                MF::IsPrimeStateActive(PrimeState::Density) &&
                MF::IsPrimeStateActive(PrimeState::VelocityY)) {
    conservatives_container[ETI(Equation::MomentumY)] =
        prime_states_container[PTI(PrimeState::VelocityY)] *
        prime_states_container[PTI(PrimeState::Density)];
  }
  if constexpr (MF::IsEquationActive(Equation::MomentumZ) &&
                MF::IsPrimeStateActive(PrimeState::Density) &&
                MF::IsPrimeStateActive(PrimeState::VelocityZ)) {
    conservatives_container[ETI(Equation::MomentumZ)] =
        prime_states_container[PTI(PrimeState::VelocityZ)] *
        prime_states_container[PTI(PrimeState::Density)];
  }
  if constexpr (MF::IsEquationActive(Equation::Gamma) &&
                MF::IsPrimeStateActive(PrimeState::gamma)) {
    conservatives_container[ETI(Equation::Gamma)] =
        GammaModelStiffenedGas::CalculateGamma(
            prime_states_container[PTI(PrimeState::gamma)]);
  }
  if constexpr (MF::IsEquationActive(Equation::Pi) &&
                MF::IsPrimeStateActive(PrimeState::pi)) {
    conservatives_container[ETI(Equation::Pi)] =
        GammaModelStiffenedGas::CalculatePi(
            prime_states_container[PTI(PrimeState::gamma)],
            prime_states_container[PTI(PrimeState::pi)]);
  }
  if constexpr (active_equations == EquationSet::GammaModel) {
    conservatives_container[ETI(Equation::Energy)] =
        GammaModelStiffenedGas::CalculateEnergy(
            prime_states_container[PTI(PrimeState::Density)],
            MF::IsPrimeStateActive(PrimeState::VelocityX)
                ? prime_states_container[PTI(PrimeState::VelocityX)]
                : 0.0,
            MF::IsPrimeStateActive(PrimeState::VelocityY)
                ? prime_states_container[PTI(PrimeState::VelocityY)]
                : 0.0,
            MF::IsPrimeStateActive(PrimeState::VelocityZ)
                ? prime_states_container[PTI(PrimeState::VelocityZ)]
                : 0.0,
            prime_states_container[PTI(PrimeState::Pressure)],
            prime_states_container[PTI(PrimeState::gamma)],
            prime_states_container[PTI(PrimeState::pi)]);
  } else {
    // The velocity equations need to be checked inside with ternary due to
    // multi-dimensionality
    if constexpr (MF::IsEquationActive(Equation::Energy) &&
                  MF::IsPrimeStateActive(PrimeState::Density) &&
                  MF::IsPrimeStateActive(PrimeState::Pressure))
      conservatives_container[ETI(Equation::Energy)] =
          eos.Energy(prime_states_container[PTI(PrimeState::Density)],
                     MF::IsPrimeStateActive(PrimeState::VelocityX)
                         ? prime_states_container[PTI(PrimeState::VelocityX)]
                         : 0.0,
                     MF::IsPrimeStateActive(PrimeState::VelocityY)
                         ? prime_states_container[PTI(PrimeState::VelocityY)]
                         : 0.0,
                     MF::IsPrimeStateActive(PrimeState::VelocityZ)
                         ? prime_states_container[PTI(PrimeState::VelocityZ)]
                         : 0.0,
                     prime_states_container[PTI(PrimeState::Pressure)]);
  }
}

template <typename ConservativesContainerType,
          typename PrimeStatesContainerType>
inline void ConservativesToPrimeStates(
    EquationOfState const &eos,
    ConservativesContainerType const &conservatives_container,
    PrimeStatesContainerType &prime_states_container) {
  if constexpr (MF::IsPrimeStateActive(PrimeState::Density) &&
                MF::IsEquationActive(Equation::Mass)) {
    prime_states_container[PTI(PrimeState::Density)] =
        conservatives_container[ETI(Equation::Mass)];
  }
  if constexpr (MF::IsEquationActive(Equation::Mass)) {
    double const one_density =
        1.0 / conservatives_container[ETI(Equation::Mass)];
    if constexpr (MF::IsPrimeStateActive(PrimeState::VelocityX) &&
                  MF::IsEquationActive(Equation::MomentumX)) {
      prime_states_container[PTI(PrimeState::VelocityX)] =
          conservatives_container[ETI(Equation::MomentumX)] * one_density;
    }
    if constexpr (MF::IsPrimeStateActive(PrimeState::VelocityY) &&
                  MF::IsEquationActive(Equation::MomentumY)) {
      prime_states_container[PTI(PrimeState::VelocityY)] =
          conservatives_container[ETI(Equation::MomentumY)] * one_density;
    }
    if constexpr (MF::IsPrimeStateActive(PrimeState::VelocityZ) &&
                  MF::IsEquationActive(Equation::MomentumZ)) {
      prime_states_container[PTI(PrimeState::VelocityZ)] =
          conservatives_container[ETI(Equation::MomentumZ)] * one_density;
    }
  }
  if constexpr (MF::IsPrimeStateActive(PrimeState::gamma) &&
                MF::IsEquationActive(Equation::Mass)) {
    prime_states_container[PTI(PrimeState::gamma)] =
        GammaModelStiffenedGas::CalculatePrimeGamma(
            conservatives_container[ETI(Equation::Gamma)]);
  }
  if constexpr (MF::IsPrimeStateActive(PrimeState::pi) &&
                MF::IsEquationActive(Equation::Mass)) {
    prime_states_container[PTI(PrimeState::pi)] =
        GammaModelStiffenedGas::CalculatePrimePi(
            prime_states_container[PTI(PrimeState::gamma)],
            conservatives_container[ETI(Equation::Pi)]);
  }
  if constexpr (active_equations == EquationSet::GammaModel) {
    prime_states_container[PTI(PrimeState::Pressure)] =
        GammaModelStiffenedGas::CalculatePressure(
            conservatives_container[ETI(Equation::Mass)],
            conservatives_container[ETI(Equation::MomentumX)],
            CC::DIM() != Dimension::One
                ? conservatives_container[ETI(Equation::MomentumY)]
                : 0.0,
            CC::DIM() == Dimension::Three
                ? conservatives_container[ETI(Equation::MomentumZ)]
                : 0.0,
            conservatives_container[ETI(Equation::Energy)],
            prime_states_container[PTI(PrimeState::gamma)],
            prime_states_container[PTI(PrimeState::pi)]);
  } else {
    // The momentum (and energy, in isentropic case) equations need to be
    // checked inside with ternary due to multi-dimensionality
    if constexpr (MF::IsPrimeStateActive(PrimeState::Pressure) &&
                  MF::IsEquationActive(Equation::Mass)) {
      prime_states_container[PTI(PrimeState::Pressure)] =
          eos.Pressure(conservatives_container[ETI(Equation::Mass)],
                       MF::IsEquationActive(Equation::MomentumX)
                           ? conservatives_container[ETI(Equation::MomentumX)]
                           : 0.0,
                       MF::IsEquationActive(Equation::MomentumY)
                           ? conservatives_container[ETI(Equation::MomentumY)]
                           : 0.0,
                       MF::IsEquationActive(Equation::MomentumZ)
                           ? conservatives_container[ETI(Equation::MomentumZ)]
                           : 0.0,
                       MF::IsEquationActive(Equation::Energy)
                           ? conservatives_container[ETI(Equation::Energy)]
                           : 0.0);
    }
    if constexpr (MF::IsPrimeStateActive(PrimeState::Temperature) &&
                  MF::IsEquationActive(Equation::Mass) &&
                  MF::IsEquationActive(Equation::Energy)) {
      prime_states_container[PTI(PrimeState::Temperature)] = eos.Temperature(
          conservatives_container[ETI(Equation::Mass)],
          MF::IsEquationActive(Equation::MomentumX)
              ? conservatives_container[ETI(Equation::MomentumX)]
              : 0.0,
          MF::IsEquationActive(Equation::MomentumY)
              ? conservatives_container[ETI(Equation::MomentumY)]
              : 0.0,
          MF::IsEquationActive(Equation::MomentumZ)
              ? conservatives_container[ETI(Equation::MomentumZ)]
              : 0.0,
          conservatives_container[ETI(Equation::Energy)]);
    }
  }
}

#endif // PRIME_STATE_HANDLER_H
