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
#ifndef FIXED_VALUE_BOUNDARY_CONDITION_H
#define FIXED_VALUE_BOUNDARY_CONDITION_H

#include "material_boundary_condition.h"
#include "boundary_constants.h"

/**
 * @brief The FixedValueBoundaryCondition class imposes a pre-defined value for each variable in all halo cells.
 */
template<BoundaryLocation LOC>
class FixedValueBoundaryCondition : public MaterialBoundaryCondition {

private:
   // fixed values for primestates and conservatives provided at the boundaries (for each material contained in the simulation)
   std::vector<std::array<double, MF::ANOE()>> const fixed_conservatives_;
   std::vector<std::array<double, MF::ANOP()>> const fixed_prime_states_;

public:
   FixedValueBoundaryCondition() = delete;
   /**
    * @brief Default constructor. See base class constructor.
    */
   explicit FixedValueBoundaryCondition( std::vector<std::array<double, MF::ANOE()>> const& fixed_conservatives,
                                         std::vector<std::array<double, MF::ANOP()>> const& fixed_prime_states ) : fixed_conservatives_( fixed_conservatives ),
                                                                                                                   fixed_prime_states_( fixed_prime_states ) {
      /** Empty besides initializer list */
   }
   ~FixedValueBoundaryCondition()                                    = default;
   FixedValueBoundaryCondition( FixedValueBoundaryCondition const& ) = delete;
   FixedValueBoundaryCondition& operator=( FixedValueBoundaryCondition const& ) = delete;
   FixedValueBoundaryCondition( FixedValueBoundaryCondition&& )                 = delete;
   FixedValueBoundaryCondition& operator=( FixedValueBoundaryCondition&& ) = delete;

   /**
    * @brief Imposes predefined values onto the respective halo cells. See base class.
    */
   void UpdateMaterialExternal( Node& node, MaterialFieldType const field_type ) const override {
      auto start_indices = BoundaryConstants<LOC>::HaloStartIndices();
      auto end_indices   = BoundaryConstants<LOC>::HaloEndIndices();

      for( auto& host_mat_block : node.GetPhases() ) {
         unsigned int const material_index( MTI( host_mat_block.first ) );
         switch( field_type ) {
            case MaterialFieldType::Conservatives: {
               for( Equation const eq : MF::ASOE() ) {
                  double( &cells )[CC::TCX()][CC::TCY()][CC::TCZ()] = host_mat_block.second.GetRightHandSideBuffer( eq );
                  for( unsigned int i = start_indices[0]; i < end_indices[0]; ++i ) {
                     for( unsigned int j = start_indices[1]; j < end_indices[1]; ++j ) {
                        for( unsigned int k = start_indices[2]; k < end_indices[2]; ++k ) {
                           cells[i][j][k] = fixed_conservatives_[material_index][ETI( eq )];
                        }
                     }
                  }
               }
            } break;
            case MaterialFieldType::PrimeStates: {
               for( PrimeState const ps : MF::ASOP() ) {
                  double( &cells )[CC::TCX()][CC::TCY()][CC::TCZ()] = host_mat_block.second.GetPrimeStateBuffer( ps );
                  for( unsigned int i = start_indices[0]; i < end_indices[0]; ++i ) {
                     for( unsigned int j = start_indices[1]; j < end_indices[1]; ++j ) {
                        for( unsigned int k = start_indices[2]; k < end_indices[2]; ++k ) {
                           cells[i][j][k] = fixed_prime_states_[material_index][PTI( ps )];
                        }
                     }
                  }
               }
            } break;
            case MaterialFieldType::Parameters: {
               throw std::runtime_error( "For the material field parameters fixed value boundary conditions are not implemented yet!" );
            } break;
            default:
               throw std::runtime_error( "Material field type not known for fixed value BC!" );
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

#endif// FIXED_VALUE_BOUNDARY_CONDITION_H
