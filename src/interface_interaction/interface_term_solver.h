/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#ifndef INTERFACE_TERM_SOLVER_H
#define INTERFACE_TERM_SOLVER_H

#include "interface_interaction/interface_term_contributions/heat_exchange_fluxes.h"
#include "interface_interaction/interface_term_contributions/interface_stress_tensor_fluxes.h"
#include "stencils/spatial_derivative_stencils/central_difference.h"
#include "stencils/spatial_derivative_stencils/fourth_order_central_difference.h"
#include "materials/material_manager.h"
#include "topology/node.h"
#include "user_specifications/numerical_setup.h"

#include "levelset/geometry/geometry_calculator_setup.h"

using GeometryCalculatorConcretization = GeometryCalculatorSetup::Concretize<geometry_calculator>::type;

/**
 * The InterfaceTermSolver class handles the solution of interface pressure, velocity, and exchange terms.
 */
class InterfaceTermSolver {

private:
   static constexpr double one_twelfth_ = 1.0 / 12.0;

   MaterialManager const& material_manager_;
   GeometryCalculatorConcretization const geometry_calculator_;
   InterfaceStressTensorFluxes const interface_stress_tensor_fluxes_;
   HeatExchangeFluxes const heat_exchange_fluxes_;

   void FillInterfaceNormalVelocityBuffer( Node const& node, double ( &u_interface_normal_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3] ) const;

   void FillDeltaApertureBuffer( Node const& node, double ( &delta_aperture_field )[CC::ICX()][CC::ICY()][CC::ICZ()][3] ) const;

public:
   InterfaceTermSolver() = delete;
   explicit InterfaceTermSolver( MaterialManager const& material_manager );
   ~InterfaceTermSolver()                            = default;
   InterfaceTermSolver( InterfaceTermSolver const& ) = delete;
   InterfaceTermSolver& operator=( InterfaceTermSolver const& ) = delete;
   InterfaceTermSolver( InterfaceTermSolver&& )                 = delete;
   InterfaceTermSolver& operator=( InterfaceTermSolver&& ) = delete;

   void SolveInterfaceInteraction( Node& node ) const;
   void WeightFaceFluxes( Node const& node, MaterialName const material,
                          double ( &face_fluxes_x )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                          double ( &face_fluxes_y )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1],
                          double ( &face_fluxes_z )[MF::ANOE()][CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1] ) const;

   void WeightVolumeForces( Node const& node, MaterialName const material,
                            double ( &volume_forces )[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()] ) const;
};

#endif// INTERFACE_TERM_SOLVER_H
