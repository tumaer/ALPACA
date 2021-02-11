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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#ifndef SYMMETRY_BOUNDARY_CONDITION_H
#define SYMMETRY_BOUNDARY_CONDITION_H

#include "material_boundary_condition.h"
#include "levelset_boundary_condition.h"
#include "boundary_constants.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief The SymmetryBoundaryCondition class implements a symmetry (mirroring) external boundary condition of the domain.
 */
template<BoundaryLocation LOC>

class SymmetryBoundaryCondition : public MaterialBoundaryCondition, public LevelsetBoundaryCondition {

   static constexpr double SymmetrySign( MaterialFieldType const field_type, unsigned int const field_index ) {
      switch( field_type ) {
         case MaterialFieldType::Conservatives:
            return SymmetrySign( MF::ASOE()[field_index] );
         case MaterialFieldType::Parameters:
            return SymmetrySign( MF::ASOPA()[field_index] );
         default:// case MaterialFieldType::PrimeStates:
            return SymmetrySign( MF::ASOP()[field_index] );
      }
   }

   static constexpr double SymmetrySign( Equation const ) { return 1.0; }
   static constexpr double SymmetrySign( PrimeState const ) { return 1.0; }
   static constexpr double SymmetrySign( Parameter const ) { return 1.0; }

   /**
    * @brief Updates the halo cells from the internal cells according to simple symmetry (same sign).
    * @param host_buffer Reference of the buffer that is to be updated.
    */
   template<class T>
   inline void UpdateSimpleSymmetry( T ( &host_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] ) const {
      auto start_indices = BoundaryConstants<LOC>::HaloStartIndices();
      auto end_indices   = BoundaryConstants<LOC>::HaloEndIndices();

      for( unsigned int i = start_indices[0]; i < end_indices[0]; ++i ) {
         for( unsigned int j = start_indices[1]; j < end_indices[1]; ++j ) {
            for( unsigned int k = start_indices[2]; k < end_indices[2]; ++k ) {
               host_buffer[i][j][k] = BoundaryConstants<LOC>::SymmetryInternalValue( host_buffer, i, j, k );
            }
         }
      }
   }

public:
   SymmetryBoundaryCondition()                                   = default;
   ~SymmetryBoundaryCondition()                                  = default;
   SymmetryBoundaryCondition( SymmetryBoundaryCondition const& ) = delete;
   SymmetryBoundaryCondition& operator=( SymmetryBoundaryCondition const& ) = delete;
   SymmetryBoundaryCondition( SymmetryBoundaryCondition&& )                 = delete;
   SymmetryBoundaryCondition& operator=( SymmetryBoundaryCondition&& ) = delete;

   /**
    * @brief Mirrors the domain values near the interface into the halo cells. See base class.
    */
   void UpdateMaterialExternal( Node& node, MaterialFieldType const field_type ) const override {
      constexpr auto start_indices = BoundaryConstants<LOC>::HaloStartIndices();
      constexpr auto end_indices   = BoundaryConstants<LOC>::HaloEndIndices();

      unsigned int const number_of_fields = MF::ANOF( field_type );
      for( auto& host_mat_block : node.GetPhases() ) {
         for( unsigned int field_index = 0; field_index < number_of_fields; ++field_index ) {
            double( &cells )[CC::TCX()][CC::TCY()][CC::TCZ()] = host_mat_block.second.GetFieldBuffer( field_type, field_index );
            for( unsigned int i = start_indices[0]; i < end_indices[0]; ++i ) {
               for( unsigned int j = start_indices[1]; j < end_indices[1]; ++j ) {
                  for( unsigned int k = start_indices[2]; k < end_indices[2]; ++k ) {
                     cells[i][j][k] = SymmetrySign( field_type, field_index ) * BoundaryConstants<LOC>::SymmetryInternalValue( cells, i, j, k );
                  }
               }
            }
         }
      }
   }

   /**
    * @brief See base class. Adjusted to symmetry condition.
    */
   void UpdateLevelsetExternal( Node& node, InterfaceBlockBufferType const buffer_type ) const override {
      /*  NH TODO this if construct should be avoided, therefore different neighbor relations
       *  in CommunicationManger needed for levelset vs. Material/Tag Halo updates.
       */
      if( node.HasLevelset() ) {
         double( &buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetBuffer( buffer_type );
         UpdateSimpleSymmetry( buffer );
      }
   }

   /**
    * @brief See base class. Adjusted to symmetry condition.
    */
   void UpdateInterfaceTagExternal( std::int8_t ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] ) const override {
      UpdateSimpleSymmetry( interface_tags );
   }

   /**
    * @brief Identifies the Location of the BoundaryCondition.
    * @return A BoundaryLocation idicating the position of the BoundaryCondition, i. e. which halo cells are updated by it.
    */
   BoundaryLocation GetLocation() const {
      return LOC;
   }
};

// Implementations of the symmetry functions for the different natural boundary locations (conservatives)
template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::East>::SymmetrySign( Equation const equation ) {
   switch( equation ) {
      case Equation::MomentumX:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::West>::SymmetrySign( Equation const equation ) {
   switch( equation ) {
      case Equation::MomentumX:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::North>::SymmetrySign( Equation const equation ) {
   switch( equation ) {
      case Equation::MomentumY:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::South>::SymmetrySign( Equation const equation ) {
   switch( equation ) {
      case Equation::MomentumY:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::Top>::SymmetrySign( Equation const equation ) {
   switch( equation ) {
      case Equation::MomentumZ:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::Bottom>::SymmetrySign( Equation const equation ) {
   switch( equation ) {
      case Equation::MomentumZ:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::East>::SymmetrySign( PrimeState const prime_state ) {
   switch( prime_state ) {
      case PrimeState::VelocityX:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::West>::SymmetrySign( PrimeState const prime_state ) {
   switch( prime_state ) {
      case PrimeState::VelocityX:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::North>::SymmetrySign( PrimeState const prime_state ) {
   switch( prime_state ) {
      case PrimeState::VelocityY:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::South>::SymmetrySign( PrimeState const prime_state ) {
   switch( prime_state ) {
      case PrimeState::VelocityY:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::Top>::SymmetrySign( PrimeState const prime_state ) {
   switch( prime_state ) {
      case PrimeState::VelocityZ:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::Bottom>::SymmetrySign( PrimeState const prime_state ) {
   switch( prime_state ) {
      case PrimeState::VelocityZ:
         return -1.0;
      default:
         return 1.0;
   }
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::East>::SymmetrySign( Parameter const ) {
   return 1.0;
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::West>::SymmetrySign( Parameter const ) {
   return 1.0;
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::North>::SymmetrySign( Parameter const ) {
   return 1.0;
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::South>::SymmetrySign( Parameter const ) {
   return 1.0;
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::Top>::SymmetrySign( Parameter const ) {
   return 1.0;
}

template<>
constexpr double SymmetryBoundaryCondition<BoundaryLocation::Bottom>::SymmetrySign( Parameter const ) {
   return 1.0;
}

#endif// SYMMETRY_BOUNDARY_CONDITION_H
