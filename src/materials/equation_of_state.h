//===----------------------- equation_of_state.h --------------------------===//
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
#ifndef EQUATION_OF_STATE_H
#define EQUATION_OF_STATE_H

#include "equation_of_state_definitions.h"
#include <vector>

/**
 * @brief The EquationOfState class defines an interface for different
 * discretizations of the equation of states to model a certain material. It
 * does not manipulate any data, but provides computational interface of type z
 * = f(x,y).
 */
class EquationOfState {

  // functions required to be implented from the derivwed class
  virtual double ComputePressure(double const mass, double const momentum_x,
                                 double const momentum_y,
                                 double const momentum_z,
                                 double const energy) const = 0;
  virtual double ComputeEnthalpy(double const mass, double const momentum_x,
                                 double const momentum_y,
                                 double const momentum_z,
                                 double const energy) const = 0;
  virtual double ComputeEnergy(double const density, double const velocity_x,
                               double const velocity_y, double const velocity_z,
                               double const pressure) const = 0;
  // Here we use the [[maybe_unused]] syntax to explicitly name the variables
  // that should be used for derived classes.
  virtual double
  ComputeTemperature([[maybe_unused]] double const mass,
                     [[maybe_unused]] double const momentum_x,
                     [[maybe_unused]] double const momentum_y,
                     [[maybe_unused]] double const momentum_z,
                     [[maybe_unused]] double const energy) const {
    return -1.0;
  }

  virtual double GetGruneisen() const = 0;
  virtual double GetGruneisen([[maybe_unused]] double const density) const {
    return GetGruneisen();
  }
  virtual double ComputePsi(double const pressure,
                            double const one_density) const = 0;
  virtual double GetGamma() const { return -1.0; }
  virtual double GetB() const { return -1.0; }
  virtual double ComputeSpeedOfSound(double const density,
                                     double const pressure) const = 0;

protected:
  // protected default constructor (can only be called from derived classes)
  explicit EquationOfState() = default;

public:
  virtual ~EquationOfState() = default;
  EquationOfState(EquationOfState const &) = delete;
  EquationOfState &operator=(EquationOfState const &) = delete;
  EquationOfState(EquationOfState &&) = delete;
  EquationOfState &operator=(EquationOfState &&) = delete;

  /**
   * @brief Computes the pressure for given input of a arbitrary mass, momentum
   * and energy according to the material equation of state.
   * @param mass The mass used for the computation.
   * @param momentum_x The momentum in x-direction used for the computation.
   * @param momentum_y The momentum in y-direction used for the computation.
   * @param momentum_z The momentum in z-direction used for the computation.
   * @param energy The energy used for the computation.
   * @return Pressure for the state imposed by the inputs of the implemented
   * material.
   */
  double Pressure(double const mass, double const momentum_x,
                  double const momentum_y, double const momentum_z,
                  double const energy) const {
    return ComputePressure(mass, momentum_x, momentum_y, momentum_z, energy);
  }

  /**
   * @brief GetEnthalpy Computes the Enthalpy based on the given inputs
   * according to the implemented equation of state and the material parameters.
   * @param mass The mass used for the computation.
   * @param momentum_x The momentum in x-direction used for the computation.
   * @param momentum_y The momentum in y-direction used for the computation.
   * @param momentum_z The momentum in z-direction used for the computation.
   * @param energy The energy used for the computation.
   * @return enthalpy for the state imposed by the inputs of the implemented
   * material.
   */
  double Enthalpy(double const mass, double const momentum_x,
                  double const momentum_y, double const momentum_z,
                  double const energy) const {
    return ComputeEnthalpy(mass, momentum_x, momentum_y, momentum_z, energy);
  }

  /**
   * @brief Computes the temperature for given input of a arbitrary mass,
   * momentum and energy according to the material equation of state.
   * @param mass The mass used for the computation.
   * @param momentum_x The momentum in x-direction used for the computation.
   * @param momentum_y The momentum in y-direction used for the computation.
   * @param momentum_z The momentum in z-direction used for the computation.
   * @param energy The energy used for the computation.
   * @return Temperature for the state imposed by the inputs of the implemented
   * material.
   * @note Returns -1.0 if the equation of state cannot supply a temperature
   * calculation.
   */
  double Temperature(double const mass, double const momentum_x,
                     double const momentum_y, double const momentum_z,
                     double const energy) const {
    return ComputeTemperature(mass, momentum_x, momentum_y, momentum_z, energy);
  }

  /**
   * @brief Computes the energy in the material for given input of density,
   * velocity and pressure.
   * @param density The density used for the computation.
   * @param velocity_x The velocity in x-direction used for the computation.
   * @param velocity_y The velocity in y-direction used for the computation.
   * @param velocity_z The velocity in z-direction used for the computation.
   * @param pressure The pressure used for the computation.
   * @return Energy for the state imposed by the inputs of the implemented
   * material.
   */
  double Energy(double const density, double const velocity_x,
                double const velocity_y, double const velocity_z,
                double const pressure) const {
    return ComputeEnergy(density, velocity_x, velocity_y, velocity_z, pressure);
  }

  /**
   * @brief Computes the speed of sound for given input of arbitrary density and
   * pressure according to the material equation of state.
   * @param density The density used for the computation.
   * @param pressure The pressure used for the computation.
   * @return Speed of sound for the state imposed by the inputs of the
   * implemented material.
   */
  double SpeedOfSound(double const density, double const pressure) const {
    return ComputeSpeedOfSound(density, pressure);
  }

  /**
   * @brief Computes the Grueneisen coefficient according to the material
   * equation of state. Dependent on material constants only so far.
   * @return Grueneisen coefficient for the implemented material.
   */
  double Gruneisen() const { return GetGruneisen(); }

  /**
   * @brief Computes the Grueneisen coefficient according to the material
   * equation of state. Dependent on material constants only so far.
   * @return Grueneisen coefficient for the implemented material.
   */
  double Gruneisen(double const density) const { return GetGruneisen(density); }

  /**
   * @brief Computes psi for given input of arbitrary density and pressure
   * according to the material equation of state in generalized form.
   * @param pressure The pressure used for the computation.
   * @param one over density ( saves costs ) .
   * @return Psi for the state imposed by the inputs of the implemented
   * material.
   */
  double Psi(double const pressure, double const one_density) const {
    return ComputePsi(pressure, one_density);
  }

  /**
   * @brief Returns the Gamma of the material.
   * @return Gamma value of the implemented material.
   */
  double Gamma() const { return GetGamma(); }

  /**
   * @brief Returns the B of the material.
   * @return B value of the implemented material.
   */
  double B() const { return GetB(); }
};

#endif // EQUATION_OF_STATE_H
