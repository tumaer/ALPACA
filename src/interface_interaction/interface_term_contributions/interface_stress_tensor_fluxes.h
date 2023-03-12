//===---------------- interface_stress_tensor_fluxes.h --------------------===//
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
#ifndef INTERFACE_STRESS_TENSOR_FLUXES_H
#define INTERFACE_STRESS_TENSOR_FLUXES_H

#include "materials/equation_of_state_definitions.h"
#include "materials/material_manager.h"
#include "topology/node.h"
#include <vector>

/**
 * @brief Class for interface stress tensor fluxes. Calculates the interface
 * fluxes caused by the stress tensor. Thus, it considers inviscid and viscous
 * contributions.
 */
class InterfaceStressTensorFluxes {

  /**
   * @brief A struct to bundle viscous flux related information of a material.
   */
  struct ViscousMaterialProperties {
    MaterialName const material_;
    double const mu_shear_;
    double const mu_bulk_;

    /**
     * @brief      Default constructor of the struct.
     *
     * @param[in]  material  The material for which viscosities are saved.
     * @param[in]  mu        A vector containing the shear and bulk viscosity of
     * the material.
     */
    ViscousMaterialProperties(MaterialName const material,
                              std::vector<double> const mu)
        : material_(material), mu_shear_(mu[0]),
          mu_bulk_(mu[1]) { /* Empty, besides initializer list */
    }
  };

  MaterialManager const &material_manager_;
  ViscousMaterialProperties const positive_material_properties_;
  ViscousMaterialProperties const negative_material_properties_;

  static constexpr double epsilon_ = std::numeric_limits<double>::epsilon();

  std::array<double, 3>
  ComputeInterfaceViscosities(double const volume_fraction) const;

  void AddInviscidPartToInterfaceStressTensor(
      Node const &node,
      double (&interface_stress_tensor_positive_material)
          [CC::ICX()][CC::ICY()][CC::ICZ()][DTI(CC::DIM())][DTI(CC::DIM())],
      double (
          &interface_stress_tensor_negative_material)[CC::ICX()][CC::ICY()]
                                                     [CC::ICZ()][DTI(CC::DIM())]
                                                     [DTI(CC::DIM())]) const;
  void AddViscousPartToInterfaceStressTensor(
      Node const &node,
      double (&interface_stress_tensor_positive_material)
          [CC::ICX()][CC::ICY()][CC::ICZ()][DTI(CC::DIM())][DTI(CC::DIM())],
      double (
          &interface_stress_tensor_negative_material)[CC::ICX()][CC::ICY()]
                                                     [CC::ICZ()][DTI(CC::DIM())]
                                                     [DTI(CC::DIM())]) const;

  void CalculateVelocityGradientAtInterface(
      Node const &node,
      double (&velocity_gradient_at_interface)[CC::ICX()][CC::ICY()][CC::ICZ()]
                                              [DTI(CC::DIM())][DTI(CC::DIM())])
      const;

  void CalculateViscousStressTensor(
      Node const &node,
      double const (&velocity_gradient)[CC::ICX()][CC::ICY()][CC::ICZ()]
                                       [DTI(CC::DIM())][DTI(CC::DIM())],
      double (&tau)[CC::ICX()][CC::ICY()][CC::ICZ()][DTI(CC::DIM())]
                   [DTI(CC::DIM())]) const;

  void AddAxisymmetricPartToViscousStressTensor(
      Node const &node, double (&tau)[CC::ICX()][CC::ICY()][CC::ICZ()]
                                     [DTI(CC::DIM())][DTI(CC::DIM())]) const;

  void AddFluxesToRightHandSide(
      Node &node,
      double const (&delta_aperture_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3],
      double const (
          &u_interface_normal_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3],
      double const (&interface_stress_tensor_positive_material)
          [CC::ICX()][CC::ICY()][CC::ICZ()][DTI(CC::DIM())][DTI(CC::DIM())],
      double const (
          &interface_stress_tensor_negative_material)[CC::ICX()][CC::ICY()]
                                                     [CC::ICZ()][DTI(CC::DIM())]
                                                     [DTI(CC::DIM())]) const;

  void ComputeRealMaterialVelocity(
      Node const &node,
      double (&real_material_velocity_x)[CC::TCX()][CC::TCY()][CC::TCZ()],
      double (&real_material_velocity_y)[CC::TCX()][CC::TCY()][CC::TCZ()],
      double (
          &real_material_velocity_z)[CC::TCX()][CC::TCY()][CC::TCZ()]) const;

public:
  InterfaceStressTensorFluxes() = delete;
  explicit InterfaceStressTensorFluxes(MaterialManager const &material_manager,
                                       MaterialName const material_positive,
                                       std::vector<double> const mu_positive,
                                       MaterialName const material_negative,
                                       std::vector<double> const mu_negative);
  ~InterfaceStressTensorFluxes() = default;
  InterfaceStressTensorFluxes(InterfaceStressTensorFluxes const &) = delete;
  InterfaceStressTensorFluxes &
  operator=(InterfaceStressTensorFluxes const &) = delete;
  InterfaceStressTensorFluxes(InterfaceStressTensorFluxes &&) = delete;
  InterfaceStressTensorFluxes &
  operator=(InterfaceStressTensorFluxes &&) = delete;

  void ComputeInterfaceFluxes(
      Node &node,
      double const (&delta_aperture_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3],
      double const (
          &u_interface_normal_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3]) const;
};

#endif // INTERFACE_STRESS_TENSOR_FLUXES_H
