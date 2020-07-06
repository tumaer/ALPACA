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
#include "source_term_solver.h"

/**
 * @brief Standard constructor using an already existing MaterialManager and the user-defined gravity.
 * @param material_manager The material manager object provides the correct equation of state for the given material.
 * @param gravity Three-dimensional array holding the gravitational pull in x-, y-, z-direction.
 */
SourceTermSolver::SourceTermSolver( MaterialManager const& material_manager, std::array<double, 3> const gravity ) :
   gravity_(gravity),
   viscous_fluxes_(material_manager),
   heat_fluxes_(material_manager),
   axisymmetric_fluxes_(),
   axisymmetric_viscous_volume_forces_( material_manager )
{
   /* Empty besides initializer list*/
}

/**
 * @brief Computes additions to the right hand side solution due to the present source terms.
 * @param mat_block The phase with its material identifier.
 * @param cell_size The cell size of the node.
 * @param node_origin_x The coordinate of the node origin in x-direction.
 * @param face_fluxes_x, face_fluxes_y, face_fluxes_z Fluxes across the cell face.
 * @param volume_forces The volume forces acting at the cell center of the node.
 */
void SourceTermSolver::Sources(std::pair<MaterialName const, Block> const& mat_block, double const cell_size, double const node_origin_x,
   double (&face_fluxes_x)[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
   double (&face_fluxes_y)[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
   double (&face_fluxes_z)[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
   double (&volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()]       ) const {

   //compute dissipative fluxes
   if constexpr( CC::ViscosityIsActive() ){
      viscous_fluxes_.ComputeFluxes( mat_block, face_fluxes_x, face_fluxes_y, face_fluxes_z, cell_size );
   }

   //compute changes due to gravity
   if constexpr( CC::GravityIsActive() ){
      gravity_.ComputeForces( mat_block.second, volume_forces );
   }

   //compute terms for axisymmetric simulations
   if constexpr( CC::Axisymmetric() ){
      axisymmetric_fluxes_.ComputeAxisymmetricContributions( mat_block.second, volume_forces, cell_size, node_origin_x );
   }

   if constexpr( CC::ViscosityIsActive() && CC::Axisymmetric() ) {
      axisymmetric_viscous_volume_forces_.ComputeForces( mat_block, volume_forces, cell_size, node_origin_x );
   }

   // Compute terms for heat exchange
   if constexpr( CC::HeatConductionActive() ) {
      heat_fluxes_.ComputeFluxes( mat_block, face_fluxes_x, face_fluxes_y, face_fluxes_z, cell_size );
   }
}
