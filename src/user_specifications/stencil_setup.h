//===------------------------ stencil_setup.h -----------------------------===//
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
#ifndef STENCIL_SETUP_H
#define STENCIL_SETUP_H

// RECONSTRUCTION_STENCIL
enum class ReconstructionStencils {
  FirstOrder,
  WENO3,
  WENOF3P,
  FourthOrderCentral,
  WENO5,
  WENO5AER,
  WENO5IS,
  WENO5Z,
  WENOAO53,
  WENO5HM,
  WENO5NU6P,
  TENO5,
  WENOCU6,
  WENO7,
  WENO9
};
constexpr ReconstructionStencils reconstruction_stencil =
    ReconstructionStencils::WENO5;
constexpr ReconstructionStencils levelset_reconstruction_stencil =
    ReconstructionStencils::WENO3;
constexpr ReconstructionStencils geometry_reconstruction_stencil =
    ReconstructionStencils::WENO3;
constexpr ReconstructionStencils viscous_fluxes_reconstruction_stencil =
    ReconstructionStencils::FourthOrderCentral;
constexpr ReconstructionStencils heat_fluxes_reconstruction_stencil =
    ReconstructionStencils::FourthOrderCentral;

// DERIVATIVE_STENCIL
enum class DerivativeStencils {
  CentralDifference,
  FourthOrderCentralDifference,
  FourthOrderCellFace,
  HOUC5
};
constexpr DerivativeStencils derivative_stencil = DerivativeStencils::HOUC5;

constexpr DerivativeStencils viscous_fluxes_derivative_stencil_cell_center =
    DerivativeStencils::FourthOrderCentralDifference;
constexpr DerivativeStencils viscous_fluxes_derivative_stencil_cell_face =
    DerivativeStencils::FourthOrderCellFace;
constexpr DerivativeStencils heat_fluxes_derivative_stencil_cell_face =
    DerivativeStencils::FourthOrderCellFace;
constexpr DerivativeStencils heat_fluxes_derivative_stencil_cell_center =
    DerivativeStencils::FourthOrderCentralDifference;

constexpr DerivativeStencils normal_calculation_derivative_stencil =
    DerivativeStencils::CentralDifference;
constexpr DerivativeStencils curvature_calculation_derivative_stencil =
    DerivativeStencils::FourthOrderCentralDifference;

#endif // STENCIL_SETUP_H
