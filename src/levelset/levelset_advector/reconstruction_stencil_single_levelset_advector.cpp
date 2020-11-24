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
#include "reconstruction_stencil_single_levelset_advector.h"

#include "enums/interface_tag_definition.h"
#include "utilities/mathematical_functions.h"
#include "enums/interface_tag_definition.h"
#include "utilities/mathematical_functions.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief Calculates the right-hand side to solve the level-set advection equation.
 *        The advection equation corresponds to equation 12 in \cite Hu2006
 * @note The implemented version differs from \cite Fedkiw2000 since the stencil average is
 *       performed only once on the gradient instead of twice on the levelset values.
 * @param node See base class.
 */
void ReconstructionStencilSingleLevelsetAdvector::AdvectImplementation( Node& node ) const {

   using ReconstructionStencil = ReconstructionStencilSetup::Concretize<levelset_reconstruction_stencil>::type;

   InterfaceBlock& interface_block = node.GetInterfaceBlock();
   double const cell_size          = node.GetCellSize();
   double const one_cell_size      = 1.0 / cell_size;

   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()]    = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
   double const( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()]               = interface_block.GetBaseBuffer( InterfaceDescription::Levelset );
   double const( &levelset_reinitialized )[CC::TCX()][CC::TCY()][CC::TCZ()] = interface_block.GetReinitializedBuffer( InterfaceDescription::Levelset );
   double( &levelset_rhs )[CC::TCX()][CC::TCY()][CC::TCZ()]                 = interface_block.GetRightHandSideBuffer( InterfaceDescription::Levelset );

   double const( &interface_velocity )[CC::TCX()][CC::TCY()][CC::TCZ()] = interface_block.GetInterfaceStateBuffer( InterfaceState::Velocity );

   std::array<double, DTI( CC::DIM() )> interface_velocity_projection;
   std::array<double, DTI( CC::DIM() )> levelset_derivative;
   for( unsigned int d = 0; d < DTI( CC::DIM() ); ++d ) {
      interface_velocity_projection[d] = 0.0;
      levelset_derivative[d]           = 0.0;
   }
   std::array<double, DTI( CC::DIM() )> increments;

   for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
            /***
             * In order to calculate the right-hand side for the level-set advection in the 2nd RK stage we need reasonable level-set values
             * in the whole extension band. Thus, it is necessary to calculate level-set advection rhs values for the extension band.
             */
            if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::ExtensionBand ) ) {

               // Calculate normal to determine x-,y- and z-component of interface_velocity
               std::array<double, 3> const normal = GetNormal( levelset_reinitialized, i, j, k );

               double const u_interface = interface_velocity[i][j][k] * one_cell_size;

               interface_velocity_projection[0] = u_interface * normal[0];
               levelset_derivative[0]           = SU::DerivativeWithUpwinding<ReconstructionStencil, Direction::X>( levelset, i, j, k, interface_velocity_projection[0], 1.0 );
               increments[0]                    = -interface_velocity_projection[0] * levelset_derivative[0];

               if constexpr( CC::DIM() != Dimension::One ) {
                  interface_velocity_projection[1] = u_interface * normal[1];
                  levelset_derivative[1]           = SU::DerivativeWithUpwinding<ReconstructionStencil, Direction::Y>( levelset, i, j, k, interface_velocity_projection[1], 1.0 );
                  increments[1]                    = -interface_velocity_projection[1] * levelset_derivative[1];
               }

               if constexpr( CC::DIM() == Dimension::Three ) {
                  interface_velocity_projection[2] = u_interface * normal[2];
                  levelset_derivative[2]           = SU::DerivativeWithUpwinding<ReconstructionStencil, Direction::Z>( levelset, i, j, k, interface_velocity_projection[2], 1.0 );
                  increments[2]                    = -interface_velocity_projection[2] * levelset_derivative[2];
               }

               levelset_rhs[i][j][k] = ConsistencyManagedSum( increments );
            }
         }//i
      }   //j
   }      //k
}
