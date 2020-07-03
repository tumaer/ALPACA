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
#ifndef STENCIL_UTILITIES_H
#define STENCIL_UTILITIES_H

#include "apply_wrapper.h"

static_assert(ReconstructionStencilSetup::Concretize<reconstruction_stencil>::type::DownstreamStencilSize() < CC::HS(), "Halo size not enough for RECONSTRUCTION_STENCIL. Increase the halo size in compile_time_constants.h!");
static_assert(DerivativeStencilSetup::Concretize<derivative_stencil>::type::DownstreamStencilSize() < CC::HS(), "Halo size not enough for DERIVATIVE_STENCIL. Increase the halo size in compile_time_constants.h!");


namespace StencilUtilities {

   template<typename S, StencilProperty P, typename T = double, Dimension DIM = CC::DIM()>
   inline T Reconstruction( std::array<double, S::StencilSize()> const& array, double const cell_size) {
      return ApplyUtilities::Apply<S, P, T>( array, cell_size );
   }

   template<typename S, StencilProperty P, Direction D, typename T = double, Dimension DIM = CC::DIM()>
   inline T Reconstruction( T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
                          , unsigned int const i, unsigned int const j, unsigned int const k
                          , double const cell_size) {
      return ApplyUtilities::Apply<S, P, D>(buffer, i, j, k, cell_size);
   }

   template<typename S, typename T = double, Dimension DIM = CC::DIM()>
   inline T ReconstructionWithUpwinding( std::array<double, S::StencilSize()> const& array, double const upwind_decision, double const cell_size) {
      if( upwind_decision >= 0 ) {
         return ApplyUtilities::Apply<S, SP::UpwindLeft>(array, cell_size);
      } else {
         return ApplyUtilities::Apply<S, SP::UpwindRight>(array, cell_size);
      }
   }

   template<typename S, Direction D, typename T = double, Dimension DIM = CC::DIM()>
   inline T ReconstructionWithUpwinding( T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
                                       , unsigned int const i, unsigned int const j, unsigned int const k
                                       , double const upwind_decision, double const cell_size) {
      if( upwind_decision >= 0 ) {
         return ApplyUtilities::Apply<S, SP::UpwindLeft, D>(buffer, i, j, k, cell_size);
      } else {
         return ApplyUtilities::Apply<S, SP::UpwindRight, D>(buffer, i, j, k, cell_size);
      }
   }


   template<typename S, typename T = double, Dimension DIM = CC::DIM()>
   inline T Derivative( std::array<double, S::StencilSize()> const& array, double const cell_size) {
      return ApplyUtilities::Apply<S, SP::Central, T>( array, cell_size );
   }

   template<typename S, Direction D, typename T = double, Dimension DIM = CC::DIM()>
   inline T Derivative( T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
                      , unsigned int const i, unsigned int const j, unsigned int const k
                      , double const cell_size) {
      switch ( S::GetStencilType() ) {
         case StencilType::Derivative: return ApplyUtilities::Apply<S, SP::Central, D>(buffer, i, j, k, cell_size);
         default: return 0.5 * ( ApplyUtilities::ApplyHjWeno<S, SP::UpwindLeft , D>( buffer, i, j, k, cell_size )
                               + ApplyUtilities::ApplyHjWeno<S, SP::UpwindRight, D>( buffer, i, j, k, cell_size ) );
      }
   }

   template<typename S, StencilProperty P, typename T = double, Dimension DIM = CC::DIM()>
   inline T Derivative( std::array<double, S::StencilSize()> const& array, double const cell_size) {
      return ApplyUtilities::Apply<S, P, T>( array, cell_size );
   }

   template<typename S, StencilProperty P, Direction D, typename T = double, Dimension DIM = CC::DIM()>
   inline T Derivative( T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
                      , unsigned int const i, unsigned int const j, unsigned int const k
                      , double const cell_size) {
      switch ( S::GetStencilType() ) {
         case StencilType::Derivative: return ApplyUtilities::Apply<S, P, D>(buffer, i, j, k, cell_size);
         default: return ApplyUtilities::ApplyHjWeno<S, P , D>( buffer, i, j, k, cell_size );
      }
   }

   template<typename S, typename T = double, Dimension DIM = CC::DIM()>
   inline T DerivativeWithUpwinding( std::array<double, S::StencilSize()> const& array
                                   , double const upwind_decision, double const cell_size) {
      if( upwind_decision >= 0 ) {
         return ApplyUtilities::Apply<S, SP::UpwindLeft, T>( array, cell_size );
      } else {
         return ApplyUtilities::Apply<S, SP::UpwindRight, T>( array, cell_size );
      }
   }

   template<typename S, Direction D, typename T = double, Dimension DIM = CC::DIM()>
   inline T DerivativeWithUpwinding( T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
                                   , unsigned int const i, unsigned int const j, unsigned int const k
                                   , double const upwind_decision, double const cell_size) {
      switch ( S::GetStencilType() ) {
         case StencilType::Derivative:
         {
            if( upwind_decision >= 0 ) {
               return ApplyUtilities::Apply<S, SP::UpwindLeft, D>(buffer, i, j, k, cell_size);
            } else {
               return ApplyUtilities::Apply<S, SP::UpwindRight, D>(buffer, i, j, k, cell_size);
            }
         }
         default:
         {
            if( upwind_decision >= 0 ) {
               return ApplyUtilities::ApplyHjWeno<S, SP::UpwindLeft , D>( buffer, i, j, k, cell_size );
            } else {
               return ApplyUtilities::ApplyHjWeno<S, SP::UpwindRight , D>( buffer, i, j, k, cell_size );
            }
         }
      }
   }

   template<typename S, typename T = double, Dimension DIM = CC::DIM()>
   inline std::array<T, 3> GradientVector(T const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                          unsigned int const i, unsigned int const j, unsigned int const k,
                                          double const cell_size) {
      return {                           Derivative<S, Direction::X, T, DIM>(buffer, i, j, k, cell_size)
             , DIM != Dimension::One ?   Derivative<S, Direction::Y, T, DIM>(buffer, i, j, k, cell_size) : 0.0
             , DIM == Dimension::Three ? Derivative<S, Direction::Z, T, DIM>(buffer, i, j, k, cell_size) : 0.0 };
   }

   template<typename S, typename T = double, Dimension DIM = CC::DIM()>
   inline std::array< std::array<T, 3>, 3> JacobianMatrix(T const (&buffer_x)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                          T const (&buffer_y)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                          T const (&buffer_z)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                          unsigned int const i, unsigned int const j, unsigned int const k,
                                                          double const cell_size) {

      return {                          GradientVector<S, T, DIM>(buffer_x, i, j, k, cell_size)
             , DIM!= Dimension::One ?   GradientVector<S, T, DIM>(buffer_y, i, j, k, cell_size) : std::array<T,3>({0.0,0.0,0.0})
             , DIM== Dimension::Three ? GradientVector<S, T, DIM>(buffer_z, i, j, k, cell_size) : std::array<T,3>({0.0,0.0,0.0})};
   }

   template<typename S, typename T = double, Dimension DIM = CC::DIM()>
   inline std::array<T, 3> Curl( T const (&buffer_x)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                 T const (&buffer_y)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                 T const (&buffer_z)[CC::TCX()][CC::TCY()][CC::TCZ()],
                                 unsigned int const i, unsigned int const j, unsigned int const k,
                                 double const cell_size ) {
      std::array< std::array<T, 3>, 3> const gradient = JacobianMatrix<S, T, DIM>(buffer_x, buffer_y, buffer_z, i, j, k, cell_size);

      return { gradient[2][1] - gradient[1][2],
               gradient[0][2] - gradient[2][0],
               gradient[1][0] - gradient[0][1] };
   }

}


namespace SU = StencilUtilities;

#endif // STENCIL_UTILITIES_H
