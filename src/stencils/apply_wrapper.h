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
#ifndef APPLY_WRAPPER_H
#define APPLY_WRAPPER_H

#include <math.h>
#include "utilities/mathematical_functions.h"

#include "user_specifications/stencil_setup.h"
#include "user_specifications/compile_time_constants.h"
#include "spatial_reconstruction_stencils/reconstruction_stencil_setup.h"
#include "spatial_derivative_stencils/derivative_stencil_setup.h"

namespace ApplyUtilities {

   /**
    * @brief Get the stencil offset and stencil sign based on whether the stencil should be evaluated UpwindLeft, UpwindRight or Central.
    * @param derivative_property Indicates whether the stencil should be evaluated UpwindLeft, UpwindRight or Central.
    * @return The stencil_offset and stencil_sign as an array.
    */
   template<StencilProperty>
   constexpr std::array<int const, 2> GetStencilParameters();
   /**
    * @brief Implementation for UpwindLeft.
    */
   template<>
   constexpr std::array<int const, 2> GetStencilParameters<SP::UpwindLeft>() {
      return {0, 1};
   }
   /**
    * @brief Implementation for UpwindRight.
    */
   template<>
   constexpr std::array<int const, 2> GetStencilParameters<SP::UpwindRight>() {
      return {1, -1};
   }
   /**
    * @brief Implementation for Central.
    */
   template<>
   constexpr std::array<int const, 2> GetStencilParameters<SP::Central>() {
      return {0, 0};
   }

   /**
    * @brief Implementation for the x-direction. See GetValueVectorFromBuffer for details.
    */
   template<typename S>
   constexpr void GetValueVectorFromBufferX( std::array<double, S::StencilSize()> & array
         , double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
         , unsigned int const i, unsigned int const j, unsigned int const k ) {
      for(unsigned int s = 0; s < S::StencilSize(); ++s) {
         array[s] = buffer[i + (-S::DownstreamStencilSize() + s)][j][k];
      }
   }
   /**
    * @brief Implementation for the y-direction. See GetValueVectorFromBuffer for details.
    */
   template<typename S>
   constexpr void GetValueVectorFromBufferY( std::array<double, S::StencilSize()> & array
         , double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
         , unsigned int const i, unsigned int const j, unsigned int const k ) {
      for(unsigned int s = 0; s < S::StencilSize(); ++s) {
         array[s] = buffer[i][j + (-S::DownstreamStencilSize() + s)][k];
      }
   }
   /**
    * @brief Implementation for the z-direction. See GetValueVectorFromBuffer for details.
    */
   template<typename S>
   constexpr void GetValueVectorFromBufferZ( std::array<double, S::StencilSize()> & array
         , double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
         , unsigned int const i, unsigned int const j, unsigned int const k ) {
      for(unsigned int s = 0; s < S::StencilSize(); ++s) {
         array[s] = buffer[i][j][k + (-S::DownstreamStencilSize() + s)];
      }
   }

   /**
    * @brief Gives the array that holds all required data from the provided buffer to apply the stencil.
    * @param array The array where the appropriate buffer data is written into (indirect return).
    * @param buffer The buffer where the array takes information from.
    * @param i,j,k The indices for which cell the stencil is applied.
    * @tparam S The used stencil.
    * @tparam D The direction in which the stencil should be evaluated.
    */
   template<typename S, Direction D>
   constexpr void GetValueVectorFromBuffer( std::array<double, S::StencilSize()> & array
         , double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
         , unsigned int const i, unsigned int const j, unsigned int const k ) {
      switch (D) {
         case Direction::X: GetValueVectorFromBufferX<S>(array, buffer, i, j, k); break;
         case Direction::Y: GetValueVectorFromBufferY<S>(array, buffer, i, j, k); break;
         default: GetValueVectorFromBufferZ<S>(array, buffer, i, j, k); break;
      }
   }

   /**
    * @brief Implementation for the x-direction. See GetDifferenceVectorFromBuffer for details.
    */
   template<typename S>
   constexpr void GetDifferenceVectorFromBufferX( std::array<double, S::StencilSize()> & array
         , double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
         , unsigned int const i, unsigned int const j, unsigned int const k
         , double const cell_size ) {
      double const one_cell_size = 1.0 / cell_size;
      for(unsigned int s = 0; s < S::StencilSize(); ++s) {
         array[s] = ( buffer[i + ( - S::DownstreamStencilSize() + s )][j][k] - buffer[i + ( -1 - S::DownstreamStencilSize() + s )][j][k] ) * one_cell_size;
      }
   }
   /**
    * @brief Implementation for the y-direction. See GetDifferenceVectorFromBuffer for details.
    */
   template<typename S>
   constexpr void GetDifferenceVectorFromBufferY( std::array<double, S::StencilSize()> & array
         , double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
         , unsigned int const i, unsigned int const j, unsigned int const k
         , double const cell_size ) {
      double const one_cell_size = 1.0 / cell_size;
      for(unsigned int s = 0; s < S::StencilSize(); ++s) {
         array[s] = ( buffer[i][j + (   - S::DownstreamStencilSize() + s)][k]
                      - buffer[i][j + (-1 - S::DownstreamStencilSize() + s)][k] ) * one_cell_size;
      }
   }
   /**
    * @brief Implementation for the z-direction. See GetDifferenceVectorFromBuffer for details.
    */
   template<typename S>
   constexpr void GetDifferenceVectorFromBufferZ( std::array<double, S::StencilSize()> & array
         , double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
         , unsigned int const i, unsigned int const j, unsigned int const k
         , double const cell_size ) {
      double const one_cell_size = 1.0 / cell_size;
      for(unsigned int s = 0; s < S::StencilSize(); ++s) {
         array[s] = ( buffer[i][j][k + (   - S::DownstreamStencilSize() + s)]
                      - buffer[i][j][k + (-1 - S::DownstreamStencilSize() + s)] ) * one_cell_size;
      }
   }

   /**
    * @brief Gives the array that holds all required data from the provided buffer to apply the stencil on difference data between two cells.
    * @param array The array where the appropriate buffer data is written into (indirect return).
    * @param buffer The buffer where the array takes information from.
    * @param i,j,k The indices for which cell the stencil is applied.
    * @tparam S The used stencil.
    * @tparam D The direction in which the stencil should be evaluated.
    */
   template<typename S, Direction D>
   constexpr void GetDifferenceVectorFromBuffer( std::array<double, S::StencilSize()> & array
         , double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()]
         , unsigned int const i, unsigned int const j, unsigned int const k
         , double const cell_size ) {
      switch (D) {
         case Direction::X: GetDifferenceVectorFromBufferX<S>(array, buffer, i, j, k, cell_size); break;
         case Direction::Y: GetDifferenceVectorFromBufferY<S>(array, buffer, i, j, k, cell_size); break;
         default: GetDifferenceVectorFromBufferZ<S>(array, buffer, i, j, k, cell_size); break;
      }
   }

   /**
    * @brief Applies the stencil on a given buffer.
    * @tparam S The used stencil.
    * @tparam P The manner in which the stencil is applied (UpwindLeft, UpwindRight or Central).
    * @tparam D The direction in which the buffer should be applied.
    * @param buffer The buffer on which the stencil is applied.
    * @param i The index in x-direction.
    * @param j The index in y-direction.
    * @param k The index in z-direction.
    * @param cell_size The cell size of the block to which the buffer belongs.
    * @return The result of the stencil.
    */
   template<typename S, StencilProperty P, Direction D>
   constexpr double Apply(double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int i, unsigned int j, unsigned int k, double const cell_size) {

      constexpr S stencil = S();

      std::array<double, S::StencilSize()> array;
      ApplyUtilities::GetValueVectorFromBuffer<S, D>(array, buffer, i, j, k);
      return stencil.template Apply<S>(array, GetStencilParameters<P>(), cell_size);
   }

   /**
    * @brief Applies the stencil on a given buffer.
    * @tparam S The used stencil.
    * @tparam P The manner in which the stencil is applied (UpwindLeft, UpwindRight or Central).
    * @tparam D The direction in which the buffer should be applied.
    * @param buffer The buffer on which the stencil is applied.
    * @param i The index in x-direction.
    * @param j The index in y-direction.
    * @param k The index in z-direction.
    * @param cell_size The cell size of the block to which the buffer belongs.
    * @return The result of the stencil.
    */
   template<typename S, StencilProperty P, Direction D>
   constexpr double ApplyHjWeno(double const (&buffer)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int i, unsigned int j, unsigned int k, double const cell_size) {

      constexpr S stencil = S();

      std::array<double, S::StencilSize()> array;
      GetDifferenceVectorFromBuffer<S, D>(array, buffer, i, j, k, cell_size);
      return stencil.template Apply<S>(array, GetStencilParameters<P>(), cell_size);
   }

   /**
    * @brief Applies the stencil on an already deduced stencil array..
    * @tparam S The used stencil.
    * @tparam P The manner in which the stencil is applied (UpwindLeft, UpwindRight or Central).
    * @tparam T The type of the output value.
    * @param array The array holding the data required for the stencil evaluation.
    * @param cell_size The cell size of the block to which the buffer belongs.
    * @return The result of the stencil.
    */
   template<typename S, StencilProperty P, typename T>
   constexpr T Apply( std::array<double, S::StencilSize()> const& array, double const cell_size) {
      constexpr S stencil = S();
      return stencil.template Apply<S>(array, ApplyUtilities::GetStencilParameters<P>(), cell_size);
   }
}

#endif // APPLY_WRAPPER_H
