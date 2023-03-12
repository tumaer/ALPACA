//===------------------------- isentropic.cpp -----------------------------===//
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
#include "isentropic.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief Constructs an isentropic-eos object with material parameters given as
 * input.
 * @param dimensional_eos_data Map containing all data for the equation of
 * state.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @note During the constructing a check is done if the required parameter
 * exists. If not an error is thrown. Furthermore, dimensionalization of each
 * value is done.
 */
Isentropic::Isentropic(
    std::unordered_map<std::string, double> const &dimensional_eos_data,
    UnitHandler const &unit_handler)
    : gamma_(GetCheckedParameter<double>(dimensional_eos_data, "gamma",
                                         "Isentropic")),
      A_(unit_handler.NonDimensionalizeValue(
          GetCheckedParameter<double>(dimensional_eos_data, "A", "Isentropic"),
          UnitType::Pressure)) {
  /* Empty besides initializer list*/
}

/**
 * @brief Computes pressure from inputs as A * rho^gamma.
 * @param mass The mass used in the computation.
 * @note Only the mass is needed in the computation in the isentropic case.
 */
double Isentropic::ComputePressure(double const mass, double const,
                                   double const, double const,
                                   double const) const {
  return A_ * std::pow(mass, gamma_);
}

/*
 * @brief Enthalpy does not matter in isentropic case.
 * @return -1.
 */
double Isentropic::ComputeEnthalpy(double const, double const, double const,
                                   double const, double const) const {
  return -1;
}

/**
 * @brief Energy does not matter in isentropic case.
 * @return -1.
 */
double Isentropic::ComputeEnergy(double const, double const, double const,
                                 double const, double const) const {
  return -1;
}

/**
 * @brief Gruneisen coefficent is not necessary for isentropic equation of
 * state.
 * @return -1.
 */
double Isentropic::GetGruneisen() const { return -1; }

/**
 * @brief Returns gamma.
 * @return gamma.
 */
double Isentropic::GetGamma() const { return gamma_; }

/**
 * @brief Background pressure is not necessary for isentropic equation of state.
 * @return -1.
 */
double Isentropic::GetB() const { return -1; }

/**
 * @brief Psi is not necessary for isentropic equation of state.
 * @param pressure .
 * @param one_density (1 / rho).
 * @return -1.
 */
double Isentropic::ComputePsi(double const, double const) const { return -1; }

/**
 * @brief Computes speed of sound from inputs as sqrt(gamma * p) / rho.
 * @param density .
 * @param pressure .
 * @return Speed of sound according to isentropic equation of state.
 */
double Isentropic::ComputeSpeedOfSound(double const density,
                                       double const pressure) const {
  return std::sqrt(gamma_ * pressure / density);
}

/**
 * @brief Provides logging information of the equation of state.
 * @param indent Number of white spaces used at the beginning of each line for
 * the logging information.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @return string with logging information.
 */
std::string Isentropic::GetLogData(unsigned int const indent,
                                   UnitHandler const &unit_handler) const {
  // string that is returned
  std::string log_string;
  // Name of the equation of state
  log_string +=
      StringOperations::Indent(indent) + "Type                 : Isentropic\n";
  // Parameters with small indentation
  log_string += StringOperations::Indent(indent) + "Gamma                : " +
                StringOperations::ToScientificNotationString(gamma_, 9) + "\n";
  log_string +=
      StringOperations::Indent(indent) + "A                    : " +
      StringOperations::ToScientificNotationString(
          unit_handler.DimensionalizeValue(A_, UnitType::Pressure), 9) +
      "\n";

  return log_string;
}
