//===-------------------------- isentropic.h ------------------------------===//
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
#ifndef ISENTROPIC_H
#define ISENTROPIC_H

#include "materials/equation_of_state.h"
#include "materials/material.h"

/**
 * @brief The Isentropic class implements a generic isentropic material, i.e.
 * material parameters must be set via caller input.
 */
class Isentropic : public EquationOfState {

  double const gamma_;
  double const A_;

  double ComputePressure(double const mass, double const momentum_x,
                         double const momentum_y, double const momentum_z,
                         double const energy) const override;
  double ComputeEnthalpy(double const mass, double const momentum_x,
                         double const momentum_y, double const momentum_z,
                         double const energy) const override;
  double ComputeEnergy(double const density, double const velocity_x,
                       double const velocity_y, double const velocity_z,
                       double const pressure) const override;
  double GetGruneisen() const override;
  double ComputePsi(double const pressure,
                    double const one_density) const override;
  double GetGamma() const override;
  double GetB() const override;
  double ComputeSpeedOfSound(double const density,
                             double const pressure) const override;

public:
  Isentropic() = delete;
  explicit Isentropic(
      std::unordered_map<std::string, double> const &dimensional_eos_data,
      UnitHandler const &unit_handler);
  virtual ~Isentropic() = default;
  Isentropic(const Isentropic &) = delete;
  Isentropic &operator=(const Isentropic &) = delete;
  Isentropic(Isentropic &&) = delete;
  Isentropic &operator=(Isentropic &&) = delete;

  std::string GetLogData(unsigned int const indent,
                         UnitHandler const &unit_handler) const;
};

#endif // ISENTROPIC_H
