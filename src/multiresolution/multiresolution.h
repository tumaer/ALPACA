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
#ifndef MULTIRESOLUTION_H
#define MULTIRESOLUTION_H

#include <cstdint>
#include "user_specifications/compile_time_constants.h"
#include "enums/norms.h"
#include "enums/remesh_identifier.h"
#include "block.h"
#include "log_writer.h"
#include "multiresolution/threshold_computer.h"

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
   ~Multiresolution() = default;
   Multiresolution( Multiresolution const& ) = delete;
   Multiresolution& operator=( Multiresolution const& ) = delete;
   Multiresolution( Multiresolution&& ) = delete;
   Multiresolution& operator=( Multiresolution&& ) = delete;

   /**
    * @brief Meta function to identify if the provided (child) node needs refinement or may be coarsened. Therefore the relative differences ("details") between the
    *        parent's prediction and the exact value of the child are compared. Details are computed according to \cite Roussel2003. In this calculation the error
    *        estimates must be adjusted by the dimensionality of the studied case. The error estimate may be computed with respect to different norms.
    * @param parent Conservative data of the parent
    * @param child Conservative data of the child.
    * @param child_id .
    * @return Remeshing decision for the provided child.
    * @tparam N The Norm used to decide whether the children should be coarsened.
    */
   template<Norm N>
   RemeshIdentifier ChildNeedsRemeshing( Block const& parent, Block const& child, std::uint64_t const child_id ) const;

   static void Average( Conservatives const& all_child_values, Conservatives& all_parent_values, std::uint64_t const child_id );
   static void AverageJumpBuffer( SurfaceBuffer const& child_values, SurfaceBuffer& parent_values,  std::uint64_t const child_id );
   static void Prediction( double const (&U_parent)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&U_child)[CC::TCX()][CC::TCY()][CC::TCZ()], const std::uint64_t child_id,
                           unsigned int const x_start = 0, unsigned int const x_count = CC::TCX(),
                           unsigned int const y_start = 0, unsigned int const y_count = CC::TCY(),
                           unsigned int const z_start = 0, unsigned int const z_count = CC::TCZ() );
   static void PropagateCutCellTagsFromChildIntoParent( std::int8_t const (&child_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                        std::int8_t (&parent_tags)[CC::TCX()][CC::TCY()][CC::TCZ()], std::uint64_t const child_id );
   static void PropagateUniformTagsFromChildIntoParent( std::int8_t const uniform_child_tag, std::int8_t (&parent_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                        std::uint64_t const child_id );
};

Multiresolution InstantiateMultiresolution( unsigned int const maximum_level, unsigned int const user_reference_level, double const user_reference_epsilon );

#endif // MULTIRESOLUTION_H
