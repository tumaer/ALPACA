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
#ifndef WALL_BOUNDARY_CONDITION_H
#define WALL_BOUNDARY_CONDITION_H

#include "material_boundary_condition.h"
#include "boundary_constants.h"

/**
 * @brief The WallBoundaryCondition class implements a wall (no-slip) external boundary condition of the domain.
 */
template<BoundaryLocation LOC>
class WallBoundaryCondition : public MaterialBoundaryCondition {

   /**
    * @brief Gives the wall sign for a specific fluid field and index.
    * @param field_type Fluid field identifier (conservatives, prime states, parameters).
    * @param field_index Index of the fluid field type.
    * @return Wall sign for the field and index.
    */
   static constexpr double WallSign( MaterialFieldType const field_type, unsigned int const field_index ) {
      switch( field_type ) {
         case MaterialFieldType::Conservatives:
            return WallSign( MF::ASOE()[field_index] );
         case MaterialFieldType::Parameters:
            return WallSign( MF::ASOPA()[field_index] );
         default: // case MaterialFieldType::PrimeStates:
            return WallSign( MF::ASOP()[field_index] );
      }
   }

   /**
    * @brief Gives the wall sign for a conservative field.
    * @param equation Conservative buffer identifier.
    * @return Wall sign.
    */
   static constexpr double WallSign( Equation const equation ) {
      switch( equation ) {
         case Equation::MomentumX:
         case Equation::MomentumY:
         case Equation::MomentumZ:
            return -1.0;
         default:
            return 1.0;
      }
   }

   /**
    * @brief Gives the wall sign for a prime state field.
    * @param prime_state Primestate buffer identifier.
    * @return Wall sign.
    */
   static constexpr double WallSign( PrimeState const prime_state ) {
      switch( prime_state ) {
         case PrimeState::VelocityX:
         case PrimeState::VelocityY:
         case PrimeState::VelocityZ:
            return -1.0;
         default:
            return 1.0;
      }
   }

   /**
    * @brief Gives the wall sign for a parameter field.
    * @param parameter Parameter buffer identifier.
    * @return Wall sign.
    *
    * @note Currently not implemented.
    */
   static constexpr double WallSign( Parameter const ) {
      // currently not needed (throw + return) due to compiler error
      throw std::runtime_error( "For the material field parameters wall value boundary conditions are not implemented yet!" );

      return 0.0;
   }

public:
   WallBoundaryCondition() = default;
   ~WallBoundaryCondition() = default;
   WallBoundaryCondition( WallBoundaryCondition const& ) = delete;
   WallBoundaryCondition& operator=( WallBoundaryCondition const& ) = delete;
   WallBoundaryCondition( WallBoundaryCondition&& ) = delete;
   WallBoundaryCondition& operator=( WallBoundaryCondition&& ) = delete;

   /**
    * @brief Mirrors the domain values near the interface into the halo cells. See base class.
    */
   void UpdateMaterialExternal( Node& node, MaterialFieldType const field_type ) const override {
      constexpr auto start_indices = BoundaryConstants<LOC>::HaloStartIndices();
      constexpr auto end_indices = BoundaryConstants<LOC>::HaloEndIndices();

      unsigned int const number_of_fields = MF::ANOF( field_type );
      for( auto& host_mat_block : node.GetPhases() ) {
         for( unsigned int field_index = 0; field_index < number_of_fields; ++field_index ) {
            double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_mat_block.second.GetFieldBuffer( field_type, field_index );
            for( unsigned int i = start_indices[0]; i < end_indices[0]; ++i ) {
               for( unsigned int j = start_indices[1]; j < end_indices[1]; ++j ) {
                  for( unsigned int k = start_indices[2]; k < end_indices[2]; ++k ) {
                     cells[i][j][k] = WallSign(field_type, field_index) * BoundaryConstants<LOC>::SymmetryInternalValue(cells, i, j, k);
                  }
               }
            }
         }
      }
   }

   /**
    * @brief Identifies the Location of the BoundaryCondition.
    * @return A BoundaryLocation indicating the position of the BoundaryCondition, i. e. which halo cells are updated by it.
    */
   BoundaryLocation GetLocation() const {
      return LOC;
   }
};

#endif // WALL_BOUNDARY_CONDITION_H
