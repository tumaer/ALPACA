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
#ifndef ZERO_GRADIENT_BOUNDARY_CONDITION_H
#define ZERO_GRADIENT_BOUNDARY_CONDITION_H

#include "boundary_constants.h"
#include "material_boundary_condition.h"
#include "levelset_boundary_condition.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief The ZeroGradientBoundaryCondition class implements a zero gradient (extending) external boundary condition of the domain.
 */
template<BoundaryLocation LOC>
class ZeroGradientBoundaryCondition : public MaterialBoundaryCondition, public LevelsetBoundaryCondition {

   /**
    * @brief Updates the halo cells from the internal cells according to the zero-gradient condition.
    * @param host_buffer Reference of the buffer that is to be updated.
    */
   template<class T>
   inline void UpdateZeroGradient( T ( &host_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] ) const {
      auto start_indices = BoundaryConstants<LOC>::HaloStartIndices();
      auto end_indices   = BoundaryConstants<LOC>::HaloEndIndices();

      for( unsigned int i = start_indices[0]; i < end_indices[0]; ++i ) {
         for( unsigned int j = start_indices[1]; j < end_indices[1]; ++j ) {
            for( unsigned int k = start_indices[2]; k < end_indices[2]; ++k ) {
               host_buffer[i][j][k] = BoundaryConstants<LOC>::ZeroGradientValue( host_buffer, i, j, k );
            }
         }
      }
   }

public:
   ZeroGradientBoundaryCondition()                                       = default;
   ~ZeroGradientBoundaryCondition()                                      = default;
   ZeroGradientBoundaryCondition( ZeroGradientBoundaryCondition const& ) = delete;
   ZeroGradientBoundaryCondition& operator=( ZeroGradientBoundaryCondition const& ) = delete;
   ZeroGradientBoundaryCondition( ZeroGradientBoundaryCondition&& )                 = delete;
   ZeroGradientBoundaryCondition&& operator=( ZeroGradientBoundaryCondition&& ) = delete;

   /**
    * @brief See base class. Imposes a zero-gradient condition at the boundary.
    */
   void UpdateMaterialExternal( Node& node, MaterialFieldType const field_type ) const override {
      unsigned int const number_of_fields = MF::ANOF( field_type );
      for( auto& host_mat_block : node.GetPhases() ) {
         for( unsigned int field_index = 0; field_index < number_of_fields; ++field_index ) {
            double( &cells )[CC::TCX()][CC::TCY()][CC::TCZ()] = host_mat_block.second.GetFieldBuffer( field_type, field_index );
            UpdateZeroGradient( cells );
         }
      }
   }

   /**
    * @brief See base class. Adjusted to zero-gradient condition.
    */
   void UpdateLevelsetExternal( Node& node, InterfaceBlockBufferType const buffer_type ) const override {
      /*  NH TODO this if construct should be avoided, therefore different neighbor relations
       *  in CommunicationManger needed for levelset vs. Material/Tag Halo updates.
       */
      if( node.HasLevelset() ) {
         double( &buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetBuffer( buffer_type );
         UpdateZeroGradient( buffer );
      }
   }

   /**
    * @brief See base class. Adjusted to zero-gradient condition.
    */
   void UpdateInterfaceTagExternal( std::int8_t ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] ) const override {
      UpdateZeroGradient( interface_tags );
   }

   /**
    * @brief Identifies the Location of the BoundaryCondition.
    * @return A BoundaryLocation indicating the position of the BoundaryCondition, i. e. which halo cells are updated by it.
    */
   BoundaryLocation GetLocation() const {
      return LOC;
   }
};

#endif// ZERO_GRADIENT_BOUNDARY_CONDITION_H
