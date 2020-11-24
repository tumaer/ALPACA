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
#include "min_iterative_levelset_reinitializer.h"
#include "utilities/buffer_operations.h"

/**
 * @brief The default constructor for a MinIterativeLevelsetReinitializer. Calls the default constructor of the base class.
 * @param halo_manager See base class.
 */
MinIterativeLevelsetReinitializer::MinIterativeLevelsetReinitializer( HaloManager& halo_manager ) : IterativeLevelsetReinitializerBase( halo_manager ) {
   // Empty Constructor, besides call of base class constructor.
}

namespace {

   /**
 * The small value used to determine whether subcell fix has to be applied or not.
 */
   constexpr double epsilon_subcell_fix_ = 1.0e-10;

   /**
 * @brief Calculate the Min-Mod value of the second derivative of the level-set field as described in \cite Min2010.
 * @param levelset The level-set field.
 * @param i The index in x-direction.
 * @param j The index in y-direction.
 * @param k The index in z-direction.
 * @param direction_0 Indicates whether we we calculate the derivative in x-direction.
 * @param direction_1 Indicates whether we we calculate the derivative in y-direction.
 * @param direction_2 Indicates whether we we calculate the derivative in z-direction.
 * @return The Min-Mod value of the second derivative.
 */
   double GetMinModOfSecondDerivative( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                       unsigned int const i, unsigned int const j, unsigned int const k,
                                       unsigned int const direction_0, unsigned int const direction_1, unsigned int const direction_2 ) {
      return MinMod( levelset[i - 1 * direction_0][j - 1 * direction_1][k - 1 * direction_2] - 2 * levelset[i][j][k] + levelset[i + 1 * direction_0][j + 1 * direction_1][k + 1 * direction_2],
                     levelset[i][j][k] - 2 * levelset[i + 1 * direction_0][j + 1 * direction_1][k + 1 * direction_2] + levelset[i + 2 * direction_0][j + 2 * direction_1][k + 2 * direction_2] );
   }

   /**
 * @brief Calculate the reinitialized level-set value for cells with a positive level-set sign for which the subcell fix applies (as described in \cite Min2010).
 * @param levelset The level-set field.
 * @param i The index in x-direction.
 * @param j The index in y-direction.
 * @param k The index in z-direction.
 * @param direction_0 Indicates whether we we calculate the derivative in x-direction.
 * @param direction_1 Indicates whether we we calculate the derivative in y-direction.
 * @param direction_2 Indicates whether we we calculate the derivative in z-direction.
 * @return The reinitialized level-set value according to the subcell fix.
 */
   double GetSubcellFixPositive( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                 unsigned int const i, unsigned int const j, unsigned int const k,
                                 unsigned int const direction_0, unsigned int const direction_1, unsigned int const direction_2 ) {

      double const levelset_xx_zero = MinMod( levelset[i - 1 * direction_0][j - 1 * direction_1][k - 1 * direction_2] - 2 * levelset[i][j][k] + levelset[i + 1 * direction_0][j + 1 * direction_1][k + 1 * direction_2],
                                              levelset[i][j][k] - 2 * levelset[i + 1 * direction_0][j + 1 * direction_1][k + 1 * direction_2] + levelset[i + 2 * direction_0][j + 2 * direction_1][k + 2 * direction_2] );

      double dx_plus = 0.0;

      double const difference = levelset[i][j][k] - levelset[i + 1 * direction_0][j + 1 * direction_1][k + 1 * direction_2];
      double const sum        = levelset[i][j][k] + levelset[i + 1 * direction_0][j + 1 * direction_1][k + 1 * direction_2];
      if( std::abs( levelset_xx_zero ) > epsilon_subcell_fix_ ) {
         double const multiplication = levelset[i][j][k] * levelset[i + 1 * direction_0][j + 1 * direction_1][k + 1 * direction_2];
         double const D              = ( 0.5 * levelset_xx_zero - sum ) * ( 0.5 * levelset_xx_zero - sum ) - 4 * multiplication;
         dx_plus                     = 0.5 + ( difference - Sign( difference ) * std::sqrt( D ) ) / ( levelset_xx_zero );
      } else {
         dx_plus = ( levelset[i][j][k] ) / ( difference );
      }

      return dx_plus;
   }

   /**
 * @brief Calculate the reinitialized level-set value for cells with a negative level-set sign for which the subcell fix applies (as described in \cite Min2010).
 * @param levelset The level-set field.
 * @param i The index in x-direction.
 * @param j The index in y-direction.
 * @param k The index in z-direction.
 * @param direction_0 Indicates whether we we calculate the derivative in x-direction.
 * @param direction_1 Indicates whether we we calculate the derivative in y-direction.
 * @param direction_2 Indicates whether we we calculate the derivative in z-direction.
 * @return The reinitialized level-set value according to the subcell fix.
 */
   double GetSubcellFixNegative( double const ( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                 unsigned int const i, unsigned int const j, unsigned int const k,
                                 unsigned int const direction_0, unsigned int const direction_1, unsigned int const direction_2 ) {

      double const levelset_xx_zero = MinMod( levelset[i - 1 * direction_0][j - 1 * direction_1][k - 1 * direction_2] - 2 * levelset[i][j][k] + levelset[i + 1 * direction_0][j + 1 * direction_1][k + 1 * direction_2],
                                              levelset[i][j][k] - 2 * levelset[i - 1 * direction_0][j - 1 * direction_1][k - 1 * direction_2] + levelset[i - 2 * direction_0][j - 2 * direction_1][k - 2 * direction_2] );

      double dx_minus = 0.0;

      double const difference = levelset[i][j][k] - levelset[i - 1 * direction_0][j - 1 * direction_1][k - 1 * direction_2];
      double const sum        = levelset[i][j][k] + levelset[i - 1 * direction_0][j - 1 * direction_1][k - 1 * direction_2];
      if( std::abs( levelset_xx_zero ) > epsilon_subcell_fix_ ) {
         double const multiplication = levelset[i][j][k] * levelset[i - 1 * direction_0][j - 1 * direction_1][k - 1 * direction_2];
         double const D              = ( 0.5 * levelset_xx_zero - sum ) * ( 0.5 * levelset_xx_zero - sum ) - 4 * multiplication;
         dx_minus                    = 0.5 + ( difference - Sign( difference ) * std::sqrt( D ) ) / ( levelset_xx_zero );
      } else {
         dx_minus = ( levelset[i][j][k] ) / ( difference );
      }

      return dx_minus;
   }

}// namespace

/**
 * @brief Reinitializes the levelset field of a node as described in \cite Min2010.
 * @param node The node with levelset block which has to be reinitialized.
 * @param levelset_type  Level-set field type that is reinitialized.
 * @return The residuum for the current node.
 */
double MinIterativeLevelsetReinitializer::ReinitializeSingleNodeImplementation( Node& node, InterfaceDescriptionBufferType const levelset_type, bool const ) const {

   std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags( levelset_type );
   InterfaceBlock& interface_block                                       = node.GetInterfaceBlock();
   double( &levelset_orig )[CC::TCX()][CC::TCY()][CC::TCZ()]             = interface_block.GetInterfaceDescriptionBuffer( levelset_type )[InterfaceDescription::Levelset];
   double const( &levelset_0_orig )[CC::TCX()][CC::TCY()][CC::TCZ()]     = interface_block.GetRightHandSideBuffer( InterfaceDescription::Levelset );

   double reinitialization_rhs[CC::TCX()][CC::TCY()][CC::TCZ()];

   double residuum = 0.0;

   double old_levelset_sign = 0.0;
   double fix_value         = 0.0;
   double increment         = 0.0;

   double derivatives[DTI( CC::DIM() )][2];
   for( unsigned int d = 0; d < DTI( CC::DIM() ); ++d ) {
      for( unsigned int e = 0; e < 2; ++e ) {
         derivatives[d][e] = 0.0;
      }
   }

   for( unsigned int i = 0; i < CC::TCX(); ++i ) {
      for( unsigned int j = 0; j < CC::TCY(); ++j ) {
         for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
            reinitialization_rhs[i][j][k] = 0.0;
         }// k
      }   // j
   }      // i

   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            old_levelset_sign = Signum( levelset_0_orig[i][j][k] );
            if( std::abs( interface_tags[i][j][k] ) > ITTI( IT::NewCutCell ) ) {

               if( levelset_orig[i][j][k] * levelset_orig[i - 1][j][k] > 0.0 || !subcell_fix_active_ ) {
                  derivatives[0][0] = levelset_orig[i][j][k] - levelset_orig[i - 1][j][k];
                  derivatives[0][0] += 0.5 * GetMinModOfSecondDerivative( levelset_orig, i - 1, j, k, 1, 0, 0 );
               } else {
                  // subcell fix
                  fix_value         = GetSubcellFixNegative( levelset_0_orig, i, j, k, 1, 0, 0 );
                  derivatives[0][0] = levelset_orig[i][j][k] / fix_value;
                  derivatives[0][0] += 0.5 * fix_value * GetMinModOfSecondDerivative( levelset_orig, i - 1, j, k, 1, 0, 0 );
               }

               if( levelset_orig[i + 1][j][k] * levelset_orig[i][j][k] > 0.0 || !subcell_fix_active_ ) {
                  derivatives[0][1] = levelset_orig[i + 1][j][k] - levelset_orig[i][j][k];
                  derivatives[0][1] -= 0.5 * GetMinModOfSecondDerivative( levelset_orig, i, j, k, 1, 0, 0 );
               } else {
                  // subcell fix
                  fix_value         = GetSubcellFixPositive( levelset_0_orig, i, j, k, 1, 0, 0 );
                  derivatives[0][1] = -levelset_orig[i][j][k] / fix_value;
                  derivatives[0][1] -= 0.5 * fix_value * GetMinModOfSecondDerivative( levelset_orig, i, j, k, 1, 0, 0 );
               }

               if constexpr( CC::DIM() != Dimension::One ) {
                  if( levelset_orig[i][j][k] * levelset_orig[i][j - 1][k] > 0.0 || !subcell_fix_active_ ) {
                     derivatives[1][0] = levelset_orig[i][j][k] - levelset_orig[i][j - 1][k];
                     derivatives[1][0] += 0.5 * GetMinModOfSecondDerivative( levelset_orig, i, j - 1, k, 0, 1, 0 );
                  } else {
                     // subcell fix
                     fix_value         = GetSubcellFixNegative( levelset_0_orig, i, j, k, 0, 1, 0 );
                     derivatives[1][0] = levelset_orig[i][j][k] / fix_value;
                     derivatives[1][0] += 0.5 * fix_value * GetMinModOfSecondDerivative( levelset_orig, i, j - 1, k, 0, 1, 0 );
                  }

                  if( levelset_orig[i][j + 1][k] * levelset_orig[i][j][k] > 0.0 || !subcell_fix_active_ ) {
                     derivatives[1][1] = levelset_orig[i][j + 1][k] - levelset_orig[i][j][k];
                     derivatives[1][1] -= 0.5 * GetMinModOfSecondDerivative( levelset_orig, i, j, k, 0, 1, 0 );
                  } else {
                     // subcell fix
                     fix_value         = GetSubcellFixPositive( levelset_0_orig, i, j, k, 0, 1, 0 );
                     derivatives[1][1] = -levelset_orig[i][j][k] / fix_value;
                     derivatives[1][1] -= 0.5 * fix_value * GetMinModOfSecondDerivative( levelset_orig, i, j, k, 0, 1, 0 );
                  }
               }

               if constexpr( CC::DIM() == Dimension::Three ) {
                  if( levelset_orig[i][j][k] * levelset_orig[i][j][k - 1] > 0.0 || !subcell_fix_active_ ) {
                     derivatives[2][0] = levelset_orig[i][j][k] - levelset_orig[i][j][k - 1];
                     derivatives[2][0] += 0.5 * GetMinModOfSecondDerivative( levelset_orig, i, j, k - 1, 0, 0, 1 );
                  } else {
                     // subcell fix
                     fix_value         = GetSubcellFixNegative( levelset_0_orig, i, j, k, 0, 0, 1 );
                     derivatives[2][0] = levelset_orig[i][j][k] / fix_value;
                     derivatives[2][0] += 0.5 * fix_value * GetMinModOfSecondDerivative( levelset_orig, i, j, k - 1, 0, 0, 1 );
                  }

                  if( levelset_orig[i][j][k + 1] * levelset_orig[i][j][k] > 0.0 || !subcell_fix_active_ ) {
                     derivatives[2][1] = levelset_orig[i][j][k + 1] - levelset_orig[i][j][k];
                     derivatives[2][1] -= 0.5 * GetMinModOfSecondDerivative( levelset_orig, i, j, k, 0, 0, 1 );
                  } else {
                     // subcell fix
                     fix_value         = GetSubcellFixPositive( levelset_0_orig, i, j, k, 0, 0, 1 );
                     derivatives[2][1] = -levelset_orig[i][j][k] / fix_value;
                     derivatives[2][1] -= 0.5 * fix_value * GetMinModOfSecondDerivative( levelset_orig, i, j, k, 0, 0, 1 );
                  }
               }

               increment = ReinitializationConstants::Dtau * old_levelset_sign * ( 1.0 - GodunovHamiltonian( derivatives, old_levelset_sign ) );
               if( ReinitializationConstants::TrackConvergence && std::abs( interface_tags[i][j][k] ) < ITTI( IT::ReinitializationBand ) ) {
                  residuum = std::max( residuum, std::abs( increment ) );
               }

               reinitialization_rhs[i][j][k] = increment;
            }
         }//k
      }   //j
   }      //i

   //add to level-set and also the residuum
   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            levelset_orig[i][j][k] += reinitialization_rhs[i][j][k];
         }//k
      }   //j
   }      //i

   return residuum;
}
