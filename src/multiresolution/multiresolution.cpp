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
#include "multiresolution/multiresolution.h"

#include <bitset>
#include "boundary_condition/boundary_specifications.h"
#include "enums/interface_tag_definition.h"
#include "utilities/string_operations.h"

/**
 * @brief Default constructor.
 * @param thresholder Thresholder used to compute the maximum allowed details in the multiresolution analysis.
 */
Multiresolution::Multiresolution( Thresholder&& thresholder ) : thresholder_( thresholder ) {
   //Empty besides initializer list
}

/**
 * @brief Averages a surface buffer into the respective parent's surface buffer.
 * @param child_buffer Buffer on the child node.
 * @param parent_buffer Buffer on the parent node.
 * @param child_id The id of the child whose values are averaged into the parent.
 * @note Adds the average to the already present value in the parent buffer.
 */
void Multiresolution::AverageJumpBuffer( SurfaceBuffer const& child_buffer, SurfaceBuffer& parent_buffer, nid_t const child_id ) {

   //compute start index of second child in parent. 0 if dimension is not considered
   static constexpr int x_index = int( CC::ICX() ) / 2;
   static constexpr int y_index = CC::DIM() != Dimension::One ? int( CC::ICY() ) / 2 : 0;
   static constexpr int z_index = CC::DIM() == Dimension::Three ? int( CC::ICZ() ) / 2 : 0;

   // clang-format off
   // std::floors ensure automated adaptation to 1D/2D
   static constexpr int first_index_look_up_table[8][6] =   {//East,West,North,South,Top, Bottom
      {     -1,       0,      -1,       0,      -1,       0}, //child_0
      {      0,      -1,      -1, y_index,      -1, z_index}, //c1
      {     -1, x_index,       0,      -1,      -1,       0}, //c2
      {x_index,      -1, y_index,      -1,      -1, z_index}, //c3
      {     -1,       0,      -1,       0,       0,      -1}, //c4
      {      0,      -1,      -1, y_index, z_index,      -1}, //c5
      {     -1, x_index,       0,      -1,       0,      -1}, //c6
      {x_index,      -1, y_index,      -1, z_index,      -1}, //c7
   };

   static constexpr int second_index_look_up_table[8][6] =  {//East,West,North,South,Top,Bottom
      {     -1,       0,      -1,       0,      -1,       0}, //child_0
      {      0,      -1,      -1,       0,      -1,       0}, //c1
      {     -1,       0,       0,      -1,      -1, z_index}, //c2
      {      0,      -1,       0,      -1,      -1, z_index}, //c3
      {     -1, x_index,      -1, y_index,       0,      -1}, //c4
      {x_index,      -1,      -1, y_index,       0,      -1}, //c5
      {     -1, x_index, y_index,      -1, z_index,      -1}, //c6
      {x_index,      -1, y_index,      -1, z_index,      -1}  //c7
   };
   // clang-format on

   for( BoundaryLocation const location : CC::ANBS() ) {
      auto const child_values = GetBoundaryJump( child_buffer, location );
      auto parent_values      = GetBoundaryJump( parent_buffer, location );

      int const index_one_start = first_index_look_up_table[PositionOfNodeAmongSiblings( child_id )][LTI( location )];
      int const index_two_start = second_index_look_up_table[PositionOfNodeAmongSiblings( child_id )][LTI( location )];

      int i_child = 0;
      int j_child = 0;

      if( index_one_start >= 0 && index_two_start >= 0 ) {
         for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
            i_child = 0;
            // std::ceil need to ensure exactly one iteration if necessary.
            for( unsigned int i = index_one_start; i < index_one_start + std::ceil( double( CC::ICY() ) / 2 ); ++i ) {
               j_child = 0;
               for( unsigned int j = index_two_start; j < index_two_start + std::ceil( double( CC::ICZ() ) / 2 ); ++j ) {
                  if constexpr( CC::DIM() == Dimension::One ) {
                     parent_values[e][i][j] += child_values[e][i_child][j_child];
                  }
                  if constexpr( CC::DIM() == Dimension::Two ) {
                     parent_values[e][i][j] += ( child_values[e][i_child][j_child] + child_values[e][i_child + 1][j_child] ) * 0.5;
                  }
                  if constexpr( CC::DIM() == Dimension::Three ) {
                     parent_values[e][i][j] += ( ( child_values[e][i_child][j_child] + child_values[e][i_child + 1][j_child + 1] ) + ( child_values[e][i_child + 1][j_child] + child_values[e][i_child][j_child + 1] ) ) * 0.25;
                  }
                  j_child += 2;
               }//j
               i_child += 2;
            }//i
         }   //eq
      }      //if
   }         //location
}

/**
 * @brief Fifth-order prediction according to \cite Kaiser2019. Applied only to cells according to the given range.
 * @param parent_values Values to be predicted into the child.
 * @param child_values Array which receives the predicted values.
 * @param child_id Id of the child that is the target of the prediction.
 * @param x_start,y_start,z_start Start index in X/Y/Z-direction. Gives the first cell number that should be filled in the child.
 * @param x_count,y_count,z_count Number of cells in X/Y/Z-direction that should be filled in the child!
 * @note No sanity checks on start or count values are done. Values in the child buffer are overriden.
 */
void Multiresolution::Prediction( double const ( &parent_values )[CC::TCX()][CC::TCY()][CC::TCZ()], double ( &child_values )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                  nid_t const child_id, unsigned int const x_start, unsigned int const x_count, unsigned int const y_start,
                                  unsigned int const y_count, unsigned int const z_start, unsigned int const z_count ) {
#ifndef PERFORMANCE
   if( x_start % 2 != 0 || y_start % 2 != 0 || z_start % 2 != 0 ) {
      throw std::logic_error( "Fatal Error in Prediction starting indices not multiple of 2" );
   }
#endif

   std::bitset<3> position( PositionOfNodeAmongSiblings( child_id ) );

   unsigned int const x_offset = position.test( 0 ) ? CC::PIOHCH() : CC::PIOLCH();
   unsigned int const y_offset = position.test( 1 ) ? CC::PIOHCH() : CC::PIOLCH();
   unsigned int const z_offset = position.test( 2 ) ? CC::PIOHCH() : CC::PIOLCH();

   // Compute the respective cells that need to be looped in the parent
   // integer divisions are automatically rounded to floor. This is intended in the expressions below.
   unsigned int const parent_x_start = ( x_start / 2 ) + x_offset;
   unsigned int const parent_y_start = CC::DIM() != Dimension::One ? ( ( y_start / 2 ) + y_offset ) : 0;
   unsigned int const parent_z_start = CC::DIM() == Dimension::Three ? ( ( z_start / 2 ) + z_offset ) : 0;
   unsigned int const parent_x_end   = parent_x_start + ( x_count / 2 );
   unsigned int const parent_y_end   = CC::DIM() != Dimension::One ? ( parent_y_start + ( y_count / 2 ) ) : 1;
   unsigned int const parent_z_end   = CC::DIM() == Dimension::Three ? ( parent_z_start + ( z_count / 2 ) ) : 1;

   unsigned int child_index_x = x_start;
   unsigned int child_index_y = CC::DIM() != Dimension::One ? y_start : 0;
   unsigned int child_index_z = CC::DIM() == Dimension::Three ? z_start : 0;

   constexpr double coefficient0 = -22.0 / 128.0;
   constexpr double coefficient1 = 3.0 / 128.0;

   constexpr double coefficient00 = coefficient0 * coefficient0;
   constexpr double coefficient01 = coefficient0 * coefficient1;
   constexpr double coefficient10 = coefficient01;
   constexpr double coefficient11 = coefficient1 * coefficient1;

   constexpr double coefficient000 = coefficient00 * coefficient0;
   constexpr double coefficient001 = coefficient00 * coefficient1;
   constexpr double coefficient010 = coefficient001;
   constexpr double coefficient011 = coefficient0 * coefficient11;
   constexpr double coefficient100 = coefficient001;
   constexpr double coefficient101 = coefficient011;
   constexpr double coefficient110 = coefficient011;
   constexpr double coefficient111 = coefficient1 * coefficient11;

   double Qx   = 0.0;
   double Qy   = 0.0;
   double Qz   = 0.0;
   double Qxy  = 0.0;
   double Qxz  = 0.0;
   double Qyz  = 0.0;
   double Qxyz = 0.0;

   // We traverse the complete parent. Loop is offsetted because stencil reaches further (fifth-order interpolation)
   /**
    * According to \cite Harten1993.
    */
   for( unsigned int i = parent_x_start; i < parent_x_end; ++i ) {
      child_index_y = y_start;
      for( unsigned int j = parent_y_start; j < parent_y_end; ++j ) {
         child_index_z = z_start;

         for( unsigned int k = parent_z_start; k < parent_z_end; ++k ) {

            // clang-format off
            //terms for 1D, 2D, 3D cases
            Qx = coefficient0 * (parent_values[i+1][j][k] - parent_values[i-1][j][k]) + coefficient1 * (parent_values[i+2][j][k] - parent_values[i-2][j][k]);

            //terms for 2D, 3D cases
            if constexpr(CC::DIM() != Dimension::One) {
               Qy = coefficient0 * (parent_values[i][j+1][k] - parent_values[i][j-1][k]) + coefficient1 * (parent_values[i][j+2][k] - parent_values[i][j-2][k]);

               Qxy = (coefficient00 * ( ( parent_values[i+1][j+1][k] + parent_values[i-1][j-1][k] ) - ( parent_values[i-1][j+1][k] + parent_values[i+1][j-1][k] ) )
                    + coefficient11 * ( ( parent_values[i+2][j+2][k] + parent_values[i-2][j-2][k] ) - ( parent_values[i-2][j+2][k] + parent_values[i+2][j-2][k] ) ))
                    +(coefficient01 * ( ( parent_values[i+1][j+2][k] + parent_values[i-1][j-2][k] ) - ( parent_values[i-1][j+2][k] + parent_values[i+1][j-2][k] ) )
                    + coefficient10 * ( ( parent_values[i+2][j+1][k] + parent_values[i-2][j-1][k] ) - ( parent_values[i-2][j+1][k] + parent_values[i+2][j-1][k] ) ));
            }

            //terms for 3D cases
            if constexpr(CC::DIM() == Dimension::Three) {
               Qz = coefficient0 * (parent_values[i][j][k+1] - parent_values[i][j][k-1]) + coefficient1 * (parent_values[i][j][k+2] - parent_values[i][j][k-2]);

               Qxz = (coefficient00 * ( ( parent_values[i+1][j][k+1] + parent_values[i-1][j][k-1]) - ( parent_values[i-1][j][k+1] +  parent_values[i+1][j][k-1]) )
                   +  coefficient11 * ( ( parent_values[i+2][j][k+2] + parent_values[i-2][j][k-2]) - ( parent_values[i-2][j][k+2] +  parent_values[i+2][j][k-2]) ))
                   + (coefficient01 * ( ( parent_values[i+1][j][k+2] + parent_values[i-1][j][k-2]) - ( parent_values[i-1][j][k+2] +  parent_values[i+1][j][k-2]) )
                   +  coefficient10 * ( ( parent_values[i+2][j][k+1] + parent_values[i-2][j][k-1]) - ( parent_values[i-2][j][k+1] +  parent_values[i+2][j][k-1]) ));

               Qyz = (coefficient00 * ( ( parent_values[i][j+1][k+1] + parent_values[i][j-1][k-1] ) - ( parent_values[i][j-1][k+1] + parent_values[i][j+1][k-1] ) )
                   +  coefficient11 * ( ( parent_values[i][j+2][k+2] + parent_values[i][j-2][k-2] ) - ( parent_values[i][j-2][k+2] + parent_values[i][j+2][k-2] ) ))
                   + (coefficient01 * ( ( parent_values[i][j+1][k+2] + parent_values[i][j-1][k-2] ) - ( parent_values[i][j-1][k+2] + parent_values[i][j+1][k-2] ) )
                   +  coefficient10 * ( ( parent_values[i][j+2][k+1] + parent_values[i][j-2][k-1] ) - ( parent_values[i][j-2][k+1] + parent_values[i][j+2][k-1] ) ));

               Qxyz = ConsistencyManagedSum(
                  coefficient000 * ConsistencyManagedSum( parent_values[i+1][j+1][k+1] - parent_values[i-1][j-1][k-1], parent_values[i-1][j-1][k+1] - parent_values[i+1][j+1][k-1],
                                                          parent_values[i+1][j-1][k-1] - parent_values[i-1][j+1][k+1], parent_values[i-1][j+1][k-1] - parent_values[i+1][j-1][k+1] )
                + coefficient111 * ConsistencyManagedSum( parent_values[i+2][j+2][k+2] - parent_values[i-2][j-2][k-2], parent_values[i-2][j-2][k+2] - parent_values[i+2][j+2][k-2],
                                                          parent_values[i+2][j-2][k-2] - parent_values[i-2][j+2][k+2], parent_values[i-2][j+2][k-2] - parent_values[i+2][j-2][k+2]  )
                , coefficient001 * ConsistencyManagedSum( parent_values[i+1][j+1][k+2] - parent_values[i-1][j-1][k-2], parent_values[i-1][j-1][k+2] - parent_values[i+1][j+1][k-2],
                                                          parent_values[i+1][j-1][k-2] - parent_values[i-1][j+1][k+2], parent_values[i-1][j+1][k-2] - parent_values[i+1][j-1][k+2]  )
                + coefficient110 * ConsistencyManagedSum( parent_values[i+2][j+2][k+1] - parent_values[i-2][j-2][k-1], parent_values[i-2][j-2][k+1] - parent_values[i+2][j+2][k-1],
                                                          parent_values[i+2][j-2][k-1] - parent_values[i-2][j+2][k+1], parent_values[i-2][j+2][k-1] - parent_values[i+2][j-2][k+1]  )
                , coefficient011 * ConsistencyManagedSum( parent_values[i+1][j+2][k+2] - parent_values[i-1][j-2][k-2], parent_values[i-1][j-2][k+2] - parent_values[i+1][j+2][k-2],
                                                          parent_values[i+1][j-2][k-2] - parent_values[i-1][j+2][k+2], parent_values[i-1][j+2][k-2] - parent_values[i+1][j-2][k+2]  )
                + coefficient100 * ConsistencyManagedSum( parent_values[i+2][j+1][k+1] - parent_values[i-2][j-1][k-1], parent_values[i-2][j-1][k+1] - parent_values[i+2][j+1][k-1],
                                                          parent_values[i+2][j-1][k-1] - parent_values[i-2][j+1][k+1], parent_values[i-2][j+1][k-1] - parent_values[i+2][j-1][k+1]  )
                , coefficient010 * ConsistencyManagedSum( parent_values[i+1][j+2][k+1] - parent_values[i-1][j-2][k-1], parent_values[i-1][j-2][k+1] - parent_values[i+1][j+2][k-1],
                                                          parent_values[i+1][j-2][k-1] - parent_values[i-1][j+2][k+1], parent_values[i-1][j+2][k-1] - parent_values[i+1][j-2][k+1]  )
                + coefficient101 * ConsistencyManagedSum( parent_values[i+2][j+1][k+2] - parent_values[i-2][j-1][k-2], parent_values[i-2][j-1][k+2] - parent_values[i+2][j+1][k-2],
                                                          parent_values[i+2][j-1][k-2] - parent_values[i-2][j+1][k+2], parent_values[i-2][j+1][k-2] - parent_values[i+2][j-1][k+2]  ) );
            }
            // We now fill the child
            //1D, 2D, 3D cases
            child_values[child_index_x  ][child_index_y  ][child_index_z  ] = parent_values[i][j][k] + DimensionAwareConsistencyManagedSum( Qx, Qy, Qz) + DimensionAwareConsistencyManagedSum( Qxy,  Qxz,  Qyz) + Qxyz;
            child_values[child_index_x+1][child_index_y  ][child_index_z  ] = parent_values[i][j][k] + DimensionAwareConsistencyManagedSum(-Qx, Qy, Qz) + DimensionAwareConsistencyManagedSum(-Qxy, -Qxz,  Qyz) - Qxyz;

            //2D, 3D cases
            if constexpr(CC::DIM() != Dimension::One) {
               child_values[child_index_x  ][child_index_y+1][child_index_z  ] = parent_values[i][j][k] + DimensionAwareConsistencyManagedSum( Qx, -Qy, Qz) + DimensionAwareConsistencyManagedSum(-Qxy,  Qxz, -Qyz) - Qxyz;
               child_values[child_index_x+1][child_index_y+1][child_index_z  ] = parent_values[i][j][k] + DimensionAwareConsistencyManagedSum(-Qx, -Qy, Qz) + DimensionAwareConsistencyManagedSum( Qxy, -Qxz, -Qyz) + Qxyz;
            }

            //3D cases
            if constexpr(CC::DIM() == Dimension::Three) {
               child_values[child_index_x  ][child_index_y  ][child_index_z+1] = parent_values[i][j][k] + ConsistencyManagedSum( Qx,  Qy, -Qz) + ConsistencyManagedSum( Qxy, -Qxz, -Qyz) - Qxyz;
               child_values[child_index_x+1][child_index_y  ][child_index_z+1] = parent_values[i][j][k] + ConsistencyManagedSum(-Qx,  Qy, -Qz) + ConsistencyManagedSum(-Qxy,  Qxz, -Qyz) + Qxyz;
               child_values[child_index_x  ][child_index_y+1][child_index_z+1] = parent_values[i][j][k] + ConsistencyManagedSum( Qx, -Qy, -Qz) + ConsistencyManagedSum(-Qxy, -Qxz,  Qyz) + Qxyz;
               child_values[child_index_x+1][child_index_y+1][child_index_z+1] = parent_values[i][j][k] + ConsistencyManagedSum(-Qx, -Qy, -Qz) + ConsistencyManagedSum( Qxy,  Qxz,  Qyz) - Qxyz;
            }
            child_index_z += 2;
            // clang-format on
         }//k-loop
         child_index_y += 2;
      }//j-loop
      child_index_x += 2;
   }//i-loop
}

/**
 * @brief Implementation of Meta function for L-infinity norm. See meta function.
 */
template<>
RemeshIdentifier Multiresolution::ChildNeedsRemeshing<Norm::Linfinity>( Block const& parent, Block const& child, nid_t const child_id ) const {

   double predicted_values[CC::TCX()][CC::TCY()][CC::TCZ()];
   double max_detail = 0.0;

   for( Equation const eq : MF::EWA() ) {
      double const( &exact_values )[CC::TCX()][CC::TCY()][CC::TCZ()] = child.GetRightHandSideBuffer( eq );
      Multiresolution::Prediction( parent.GetRightHandSideBuffer( eq ), predicted_values, child_id );
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               max_detail = std::max( max_detail, std::abs( exact_values[i][j][k] - predicted_values[i][j][k] ) / std::abs( exact_values[i][j][k] ) );
            }// Loop : k
         }   // Loop : j
      }      // Loop : i
   }         // Loop : Equation

   return RemeshingDecision( max_detail, LevelOfNode( child_id ) );
}

/**
 * @brief Implementation of meta function for L-one norm in three dimensions. See meta function.
 */
template<>
RemeshIdentifier Multiresolution::ChildNeedsRemeshing<Norm::Lone>( Block const& parent, Block const& child, nid_t const child_id ) const {

   double predicted_values[CC::TCX()][CC::TCY()][CC::TCZ()];
   double max_detail = 0.0;
   double error_norm = 0.0;

   double const one_number_of_cells = 1.0 / ( CC::TCX() * CC::TCY() * CC::TCZ() );

   for( Equation const eq : MF::EWA() ) {
      double const( &exact_values )[CC::TCX()][CC::TCY()][CC::TCZ()] = child.GetRightHandSideBuffer( eq );
      Multiresolution::Prediction( parent.GetRightHandSideBuffer( eq ), predicted_values, child_id );
      error_norm = 0.0;
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               error_norm += std::abs( exact_values[i][j][k] - predicted_values[i][j][k] ) / std::abs( exact_values[i][j][k] );
            }
         }
      }
      max_detail = std::max( max_detail, error_norm * one_number_of_cells );// We save the maximum values over the eq-loops.
   }                                                                        // Loop : eq

   return RemeshingDecision( max_detail, LevelOfNode( child_id ) );
}

/**
 * @brief Implementation of meta function for L-two norm in three dimensions. See meta function.
 */
template<>
RemeshIdentifier Multiresolution::ChildNeedsRemeshing<Norm::Ltwo>( Block const& parent, Block const& child, nid_t const child_id ) const {

   double predicted_values[CC::TCX()][CC::TCY()][CC::TCZ()];
   double max_detail = 0.0;
   double error_norm = 0.0;

   double const one_number_of_cells = 1.0 / ( CC::TCX() * CC::TCY() * CC::TCZ() );

   for( Equation const eq : MF::EWA() ) {
      double const( &exact_values )[CC::TCX()][CC::TCY()][CC::TCZ()] = child.GetRightHandSideBuffer( eq );
      Multiresolution::Prediction( parent.GetRightHandSideBuffer( eq ), predicted_values, child_id );
      error_norm = 0.0;
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               double const temp = ( exact_values[i][j][k] - predicted_values[i][j][k] ) / exact_values[i][j][k];
               error_norm        = error_norm + temp * temp;
            }
         }
      }
      max_detail = std::max( max_detail, std::sqrt( error_norm ) * one_number_of_cells );
   }// loop : eq

   return RemeshingDecision( max_detail, LevelOfNode( child_id ) );
}

/**
 * @brief Gives the remeshing decision for the given (relevant) detail and the level of the analysis.
 * @param detail The relevant detail according to the used error norm.
 * @param level The level on which the decision is made.
 * @return Remeshing decision.
 */
RemeshIdentifier Multiresolution::RemeshingDecision( double const detail, unsigned int const level ) const {

   /* The 32 pops out of \cite Harten 1995. = 2^(p+1). With p being the order of the prediction minus 1. Our prediction is of fifth order.
    * If the prediction order is changed this function needs to be changed. As this seems an unlikely development, we keep it hard-coded.
    */
   constexpr double fifth_order_coefficient = 32.0;

   double const epsilon_coarsen    = thresholder_.ThresholdOnLevel( level );
   double const epsilon_refinement = epsilon_coarsen * fifth_order_coefficient;

   if( detail <= epsilon_coarsen ) {
      return RemeshIdentifier::Coarse;
   } else if( detail >= epsilon_refinement ) {
      return RemeshIdentifier::Refine;
   } else {
      return RemeshIdentifier::Neutral;
   }
}

/**
 * @brief Averageing-operator equivalent for interface tags.
 * @param child_tags The child's interface tag buffer.
 * @param parent_tags The parent's interface tags buffer to receive the "averaged" values.
 * @param child_id The id of the child.
 * @note Not a real averaging. Only importance is that the cut cell markers are correctly propagated. Therefore a multiplication is (mis-)used.
 */
void Multiresolution::PropagateCutCellTagsFromChildIntoParent( std::int8_t const ( &child_tags )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                               std::int8_t ( &parent_tags )[CC::TCX()][CC::TCY()][CC::TCZ()], nid_t const child_id ) {

   unsigned int const x_start = xyz_look_up_table_[0][PositionOfNodeAmongSiblings( child_id )];
   unsigned int const y_start = CC::DIM() != Dimension::One ? xyz_look_up_table_[1][PositionOfNodeAmongSiblings( child_id )] : 0;
   unsigned int const z_start = CC::DIM() == Dimension::Three ? xyz_look_up_table_[2][PositionOfNodeAmongSiblings( child_id )] : 0;

   unsigned int const x_end = x_start + CC::PSOCICX();
   unsigned int const y_end = y_start + CC::PSOCICY();
   unsigned int const z_end = z_start + CC::PSOCICZ();

   constexpr unsigned int i_child_start = CC::FICX();
   constexpr unsigned int j_child_start = CC::DIM() != Dimension::One ? CC::FICY() : 0;
   constexpr unsigned int k_child_start = CC::DIM() == Dimension::Three ? CC::FICZ() : 0;

   unsigned int i_child = i_child_start;
   unsigned int j_child = j_child_start;
   unsigned int k_child = k_child_start;

   int overflow_safe_buffer = 0.0;

   for( unsigned int i = x_start; i < x_end; ++i ) {
      j_child = j_child_start;
      for( unsigned int j = y_start; j < y_end; ++j ) {
         k_child = k_child_start;
         for( unsigned int k = z_start; k < z_end; ++k ) {
            //ATTENTION: in the following, we multiply bools whith each other. This is fine as they are implicitely casted to int.
            overflow_safe_buffer = ( std::abs( child_tags[i_child][j_child][k_child] ) > ITTI( IT::NewCutCell ) ) * ( std::abs( child_tags[i_child + 1][j_child][k_child] ) > ITTI( IT::NewCutCell ) )
#if DIMENSION == 1
#elif DIMENSION == 2
                                   * ( std::abs( child_tags[i_child][j_child + 1][k_child] ) > ITTI( IT::NewCutCell ) ) * ( std::abs( child_tags[i_child + 1][j_child + 1][k_child] ) > ITTI( IT::NewCutCell ) )
#else
                                   * ( std::abs( child_tags[i_child][j_child + 1][k_child] ) > ITTI( IT::NewCutCell ) ) * ( std::abs( child_tags[i_child + 1][j_child + 1][k_child] ) > ITTI( IT::NewCutCell ) ) * ( std::abs( child_tags[i_child][j_child][k_child + 1] ) > ITTI( IT::NewCutCell ) ) * ( std::abs( child_tags[i_child + 1][j_child][k_child + 1] ) > ITTI( IT::NewCutCell ) ) * ( std::abs( child_tags[i_child][j_child + 1][k_child + 1] ) > ITTI( IT::NewCutCell ) ) * ( std::abs( child_tags[i_child + 1][j_child + 1][k_child + 1] ) > ITTI( IT::NewCutCell ) )
#endif
                  ;
            // A cut cell is shown here as a ZERO, as we compare it to be larger than 1.
            // Mixed signs in the considered child cells are only possible if there is also a ZERO, hence overflow_safe_buffer would be ZERO as well.
            // It is therefore safe to use both the sign of overflow_safe_buffer and an arbitrary child cell to determine the sign (or ZERO) of the parent value.
            parent_tags[i][j][k] = static_cast<std::int8_t>( Signum( overflow_safe_buffer ) ) * Signum( child_tags[i_child][j_child][k_child] ) * ITTI( IT::BulkPhase );
            k_child += 2;
         }
         j_child += 2;
      }
      i_child += 2;
   }
}

/**
 * @brief Populates the parent's interface tags with the given uniform child value.
 * @param child_tag The uniform interface tag.
 * @param parent_tags The parent's tag array to be updated based on the child's tag.
 * @param child_id The id of the child.
 */
void Multiresolution::PropagateUniformTagsFromChildIntoParent( std::int8_t const child_tag, std::int8_t ( &parent_tags )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                               nid_t const child_id ) {

   unsigned int const x_start = xyz_look_up_table_[0][PositionOfNodeAmongSiblings( child_id )];
   unsigned int const y_start = CC::DIM() != Dimension::One ? xyz_look_up_table_[1][PositionOfNodeAmongSiblings( child_id )] : 0;
   unsigned int const z_start = CC::DIM() == Dimension::Three ? xyz_look_up_table_[2][PositionOfNodeAmongSiblings( child_id )] : 0;

   unsigned int const x_end = x_start + CC::PSOCICX();
   unsigned int const y_end = y_start + CC::PSOCICY();
   unsigned int const z_end = z_start + CC::PSOCICZ();

   for( unsigned int i = x_start; i < x_end; ++i ) {
      for( unsigned int j = y_start; j < y_end; ++j ) {
         for( unsigned int k = z_start; k < z_end; ++k ) {
            parent_tags[i][j][k] = child_tag;
         }
      }
   }
}
