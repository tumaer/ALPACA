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
#include "space_solver.h"

#include <cmath>
#include "utilities/mathematical_functions.h"

/**
 * @brief Standard constructor using an already existing MaterialManager and the user-defined gravity.
 * @param material_manager The material manager object provides the correct equation of state for the given material.
 * @param gravity Three-dimensional array holding the gravitational pull in x-, y-, z-direction.
 */
SpaceSolver::SpaceSolver( MaterialManager const& material_manager, std::array<double, 3> const gravity ) :
   eigendecomposition_calculator_( material_manager ),
   riemann_solver_( material_manager,eigendecomposition_calculator_ ),
   source_term_solver_( material_manager, gravity ),
   interface_term_solver_( material_manager ),
   levelset_advector_()
{
   /* Empty besides initializer list*/
}

/**
 * @brief Computes right side of the underlying system of equations (including source terms).
 * @param node The node under consideration.
 */
void SpaceSolver::UpdateFluxes( Node& node ) const {

   double face_fluxes_x[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];
   double face_fluxes_y[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];
   double face_fluxes_z[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];

   double volume_forces[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()];

   for( auto& phase : node.GetPhases() ) {
      for( const Equation eq : MF::ASOE() ) {
         double (&rhs_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.second.GetRightHandSideBuffer( eq );
         for( unsigned int i = 0; i < CC::TCX(); ++i ) {
            for( unsigned int j = 0; j < CC::TCY(); ++j ) {
               for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                  rhs_buffer[i][j][k] = 0.0;
               } //k
            } //j
         } //i
      } //equation
   } //phases

   // For multi-phase + Lmax nodes (which have a levelset)
   if( node.HasLevelset() ) {

      double (&levelset_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetRightHandSideBuffer(InterfaceDescription::Levelset);
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               levelset_rhs[i][j][k] = 0.0;
            } //k
         } //j
      } //i

      // Solve interface Riemann problem to obtain interface velocity and interface exchange terms
      interface_term_solver_.SolveInterfaceInteraction( node );
      // compute levelset rhs (from levelset and interface velocity)
      levelset_advector_.Advect( node );
   }

   double const cell_size = node.GetCellSize();
   double const one_cell_size = 1.0 / cell_size;

   for( auto& phase : node.GetPhases() ) {
      // RHS-buffers have to be reset for each phase!
      for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
         for( unsigned int i = 0; i < CC::ICX()+1; ++i ) {
            for( unsigned int j = 0; j < CC::ICY()+1; ++j ) {
               for( unsigned int k = 0; k < CC::ICZ()+1; ++k ) {
                  face_fluxes_x[e][i][j][k] = 0.0;
                  face_fluxes_y[e][i][j][k] = 0.0;
                  face_fluxes_z[e][i][j][k] = 0.0;
               } //k
            } //j
         } //i
      } //equation

      for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
         for( unsigned int i = 0; i < CC::ICX(); ++i ) {
            for( unsigned int j = 0; j < CC::ICY(); ++j ) {
               for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                  volume_forces[e][i][j][k] = 0.0;
               } //k
            } //j
         } //i
      } //equation

      // Determine cell face fluxes unsing a Riemann solver
      if constexpr( CC::InviscidExchangeActive() ) {
         riemann_solver_.Update( phase, cell_size, face_fluxes_x, face_fluxes_y, face_fluxes_z );
      }

      // Determine source terms
      source_term_solver_.Sources( phase, cell_size, node.GetBlockCoordinateX(), face_fluxes_x, face_fluxes_y, face_fluxes_z, volume_forces );

      if( node.HasLevelset() ) {
         interface_term_solver_.WeightFaceFluxes( node, phase.first, face_fluxes_x, face_fluxes_y, face_fluxes_z );
         interface_term_solver_.WeightVolumeForces( node, phase.first, volume_forces );
      }

      //update cells due to fluxes
      for( const Equation eq : MF::ASOE() ) {
         auto const e = ETI(eq);
         double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = phase.second.GetRightHandSideBuffer( eq );
         for( unsigned int i = 0; i < CC::ICX(); ++i ) {
            for( unsigned int j = 0; j < CC::ICY(); ++j ) {
               for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
                  // We need to add here because the RHS buffer might already be filled in cut cells for single-level-set simulations
                  cells[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()]  += volume_forces[e][i][j][k]
                     + DimensionAwareConsistencyManagedSum( face_fluxes_x[e][i]  [j+1][k+1] - face_fluxes_x[e][i+1][j+1][k+1],
                                                            face_fluxes_y[e][i+1][j]  [k+1] - face_fluxes_y[e][i+1][j+1][k+1],
                                                            face_fluxes_z[e][i+1][j+1][k]   - face_fluxes_z[e][i+1][j+1][k+1]
                                                          ) * one_cell_size;
               } //k
            } //j
         } //i
      } //equation

      //save boundary fluxes for correction at jump boundary conditions
      double   (&boundary_fluxes_west)[MF::ANOE()][CC::ICY()][CC::ICZ()] = phase.second.GetBoundaryJumpFluxes( BoundaryLocation::West );
      double   (&boundary_fluxes_east)[MF::ANOE()][CC::ICY()][CC::ICZ()] = phase.second.GetBoundaryJumpFluxes( BoundaryLocation::East );
      double  (&boundary_fluxes_south)[MF::ANOE()][CC::ICY()][CC::ICZ()] = phase.second.GetBoundaryJumpFluxes( BoundaryLocation::South );
      double  (&boundary_fluxes_north)[MF::ANOE()][CC::ICY()][CC::ICZ()] = phase.second.GetBoundaryJumpFluxes( BoundaryLocation::North );
      double (&boundary_fluxes_bottom)[MF::ANOE()][CC::ICY()][CC::ICZ()] = phase.second.GetBoundaryJumpFluxes( BoundaryLocation::Bottom );
      double    (&boundary_fluxes_top)[MF::ANOE()][CC::ICY()][CC::ICZ()] = phase.second.GetBoundaryJumpFluxes( BoundaryLocation::Top );
      for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
         for( unsigned int i = 0; i < CC::ICY(); ++i ) {
            for( unsigned int j = 0; j < CC::ICZ(); ++j ) {
               boundary_fluxes_west[e][i][j]   += face_fluxes_x[e][0][i+1][j+1];
               boundary_fluxes_east[e][i][j]   += face_fluxes_x[e][CC::ICX()][i+1][j+1];
               boundary_fluxes_south[e][i][j]  += face_fluxes_y[e][i+1][0][j+1];
               boundary_fluxes_north[e][i][j]  += face_fluxes_y[e][i+1][CC::ICY()][j+1];
               boundary_fluxes_bottom[e][i][j] += face_fluxes_z[e][i+1][j+1][0];
               boundary_fluxes_top[e][i][j]    += face_fluxes_z[e][i+1][j+1][CC::ICZ()];
            }
         }
      }
   } //phases
}

/**
 * @brief Computes the maximum eigenvalue in the given node.
 * @param block The block for which the eigenvalues are to be computed.
 * @param eigenvalues Indirect return parameter.
 */
void SpaceSolver::ComputeMaxEigenvaluesForPhase( std::pair<const MaterialName, Block> const& mat_block, double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()] ) const {
   eigendecomposition_calculator_.ComputeMaxEigenvaluesOnBlock( mat_block, eigenvalues );
}

/**
 * @brief Stores the given (GLF) eigenvalues for later usage in a globally valid location.
 * @param eigenvalues Values to be stored.
 */
void SpaceSolver::SetFluxFunctionGlobalEigenvalues( double (&eigenvalues)[DTI(CC::DIM())][MF::ANOE()] ) const {
   eigendecomposition_calculator_.SetGlobalEigenvalues( eigenvalues );
}
