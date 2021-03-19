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
#ifndef MULTIRESOLUTION_H
#define MULTIRESOLUTION_H

#include "user_specifications/compile_time_constants.h"
#include "enums/norms.h"
#include "enums/remesh_identifier.h"
#include "block_definitions/block.h"
#include "input_output/log_writer/log_writer.h"
#include "multiresolution/threshold_computer.h"
#include "topology/id_information.h"
#include "utilities/mathematical_functions.h"

namespace {
   constexpr unsigned int xyz_look_up_table_[3][8] = {
         { CC::PIOLCFICX(), CC::PIOHCFICX(), CC::PIOLCFICX(), CC::PIOHCFICX(), CC::PIOLCFICX(), CC::PIOHCFICX(), CC::PIOLCFICX(), CC::PIOHCFICX() },// X-Indecies
         { CC::PIOLCFICY(), CC::PIOLCFICY(), CC::PIOHCFICY(), CC::PIOHCFICY(), CC::PIOLCFICY(), CC::PIOLCFICY(), CC::PIOHCFICY(), CC::PIOHCFICY() },// Y-Indecies
         { CC::PIOLCFICZ(), CC::PIOLCFICZ(), CC::PIOLCFICZ(), CC::PIOLCFICZ(), CC::PIOHCFICZ(), CC::PIOHCFICZ(), CC::PIOHCFICZ(), CC::PIOHCFICZ() } };
}

/**
 * @brief The Multiresolution class provides prediction and averaging methods for field and interface tags as well as more high-level helper functions based
 * on either of the two.
 */
class Multiresolution {

   Thresholder const thresholder_;

   RemeshIdentifier RemeshingDecision( double const detail, unsigned int const level ) const;

public:
   Multiresolution() = delete;
   explicit Multiresolution( Thresholder&& thresholder );
   ~Multiresolution()                        = default;
   Multiresolution( Multiresolution const& ) = delete;
   Multiresolution& operator=( Multiresolution const& ) = delete;
   Multiresolution( Multiresolution&& )                 = delete;
   Multiresolution& operator=( Multiresolution&& ) = delete;

   /**
    * @brief Meta function to identify if the provided (child) node needs refinement or may be coarsened. Therefore the relative differences ("details") between the
    *        parent's prediction and the exact value of the child are compared. Details are computed according to \cite Roussel2003. In this calculation the error
    *        estimates must be adjusted by the dimensionality of the studied case. The error estimate may be computed with respect to different norms.
    * @param parent Conservative data of the parent.
    * @param child Conservative data of the child.
    * @param child_id The id of the child node.
    * @return Remeshing decision for the provided child.
    * @tparam N The Norm used to decide whether the children should be coarsened.
    */
   template<Norm N>
   RemeshIdentifier ChildNeedsRemeshing( Block const& parent, Block const& child, nid_t const child_id ) const;

   /**
    * @brief Averages the child values into the parent, i.e. conservative average of the eight (in 3D) child cells that make up one parent cell.
    * @param child_buffer The child's buffer.
    * @param parent_buffer The parent's buffer to receive the averaged values.
    * @param child_id The id of the child.
    * @note Overrides the values in the parent buffer.
    */
   template<typename BufferType>
   static void Average( BufferType const& child_buffer, BufferType& parent_buffer, nid_t const child_id ) {

      unsigned int const x_start = xyz_look_up_table_[0][PositionOfNodeAmongSiblings( child_id )];
      unsigned int const y_start = CC::DIM() != Dimension::One ? xyz_look_up_table_[1][PositionOfNodeAmongSiblings( child_id )] : 0;
      unsigned int const z_start = CC::DIM() == Dimension::Three ? xyz_look_up_table_[2][PositionOfNodeAmongSiblings( child_id )] : 0;

      unsigned int const x_end = x_start + CC::PSOCICX();
      unsigned int const y_end = y_start + CC::PSOCICY();
      unsigned int const z_end = z_start + CC::PSOCICZ();

      unsigned int const i_child_start = CC::FICX();
      unsigned int const j_child_start = CC::DIM() != Dimension::One ? CC::FICY() : 0;
      unsigned int const k_child_start = CC::DIM() == Dimension::Three ? CC::FICZ() : 0;

      unsigned int i_child = i_child_start;
      unsigned int j_child = j_child_start;
      unsigned int k_child = k_child_start;

      for( size_t field_index = 0; field_index < BufferType::GetNumberOfFields(); ++field_index ) {
         auto const& child_values = child_buffer[field_index];
         auto& parent_values      = parent_buffer[field_index];

         i_child = i_child_start;
         for( unsigned int i = x_start; i < x_end; ++i ) {
            j_child = j_child_start;
            for( unsigned int j = y_start; j < y_end; ++j ) {
               k_child = k_child_start;
               for( unsigned int k = z_start; k < z_end; ++k ) {
                  if constexpr( CC::DIM() == Dimension::One ) {
                     parent_values[i][j][k] = ( child_values[i_child][j_child][k_child] + child_values[i_child + 1][j_child][k_child] ) * 0.5;
                  }
                  if constexpr( CC::DIM() == Dimension::Two ) {
                     parent_values[i][j][k] = ( ( child_values[i_child][j_child][k_child] + child_values[i_child + 1][j_child + 1][k_child] ) + ( child_values[i_child + 1][j_child][k_child] + child_values[i_child][j_child + 1][k_child] ) ) * 0.25;
                  }
                  if constexpr( CC::DIM() == Dimension::Three ) {
                     parent_values[i][j][k] = 0.125 * ConsistencyManagedSum(
                                                            child_values[i_child][j_child][k_child] + child_values[i_child + 1][j_child + 1][k_child + 1],
                                                            child_values[i_child][j_child][k_child + 1] + child_values[i_child + 1][j_child + 1][k_child],
                                                            child_values[i_child][j_child + 1][k_child] + child_values[i_child + 1][j_child][k_child + 1],
                                                            child_values[i_child][j_child + 1][k_child + 1] + child_values[i_child + 1][j_child][k_child] );
                  }
                  k_child += 2;
               }//k
               j_child += 2;
            }//j
            i_child += 2;
         }//i
      }   //eq
   }

   static void AverageJumpBuffer( SurfaceBuffer const& child_values, SurfaceBuffer& parent_values, nid_t const child_id );
   static void Prediction( double const ( &U_parent )[CC::TCX()][CC::TCY()][CC::TCZ()], double ( &U_child )[CC::TCX()][CC::TCY()][CC::TCZ()], nid_t const child_id,
                           unsigned int const x_start = 0, unsigned int const x_count = CC::TCX(),
                           unsigned int const y_start = 0, unsigned int const y_count = CC::TCY(),
                           unsigned int const z_start = 0, unsigned int const z_count = CC::TCZ() );
   static void PropagateCutCellTagsFromChildIntoParent( std::int8_t const ( &child_tags )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                        std::int8_t ( &parent_tags )[CC::TCX()][CC::TCY()][CC::TCZ()], nid_t const child_id );
   static void PropagateUniformTagsFromChildIntoParent( std::int8_t const uniform_child_tag, std::int8_t ( &parent_tags )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                        nid_t const child_id );
};

#endif// MULTIRESOLUTION_H
