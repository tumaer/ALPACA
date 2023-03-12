//===-------------------------- field_enums.h -----------------------------===//
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
#ifndef FIELD_ENUMS_H
#define FIELD_ENUMS_H

/**
 * @brief Unique identifier for all possible conservative equations in arbitrary
 * order.
 */
enum class EquationPool {
  Mass,
  MomentumX,
  MomentumY,
  MomentumZ,
  Energy,
  Gamma,
  Pi, // Example
};

/**
 * @brief Unique Identifier for all possible prime states in arbitrary order.
 */
enum class PrimeStatePool {
  Density,
  Pressure,
  VelocityX,
  VelocityY,
  VelocityZ,
  Temperature,
  gamma,
  pi // Example
};
/**
 * @brief Unique Identifier for all possible interface quantities in arbitrary
 * order.
 */
enum class InterfaceStatePool { Velocity, PressurePositive, PressureNegative };

/**
 * @brief Unique Identifier for all possible interface descriptions variables in
 * arbitrary order.
 * @note  Every member has to be added to the InterfaceDescription enumeration
 * as well.
 */
enum class InterfaceDescriptionPool { Levelset, VolumeFraction };

/**
 * @brief Unique Identifier for all possible interface parameters in arbitrary
 * order.
 * @note  Every member has to be added to the InterfaceParameter enumeration as
 * well.
 */
enum class InterfaceParameterPool { SurfaceTensionCoefficient };

/**
 * @brief Unique Identifier for all possible parameters in arbitrary order.
 * @note  Every member has to be added to the Parameters enumeration as well.
 */
enum class ParameterPool {
  // viscosity parameters
  ShearViscosity,
  // conductivity parameters
  ThermalConductivity
};

#endif // FIELD_ENUMS_H
