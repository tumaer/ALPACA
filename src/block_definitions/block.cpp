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
#include "block_definitions/block.h"
#include <stdexcept>
#include "utilities/buffer_operations.h"

/**
 * @brief Standard constructor, creates a Block of the provided material. Initializes all buffers
 *        to zero (important with first touch rule on distributed Memory Machines)!
 */
Block::Block() {

   // We initialize a Block with zero in all its buffers
   BO::SetFieldBuffer( GetAverageBuffer(), 0.0 );
   BO::SetFieldBuffer( GetRightHandSideBuffer(), 0.0 );
   BO::SetFieldBuffer( GetInitialBuffer(), 0.0 );
   BO::SetFieldBuffer( GetPrimeStateBuffer(), 0.0 );

   // Only reset buffer of parameter if they are present
   if constexpr( CC::ParameterModelActive() ) {
      BO::SetFieldBuffer( GetParameterBuffer(), 0.0 );
   }

   for( BoundaryLocation const location : CC::NBS() ) {
      ResetJumpFluxes( location );
   }

   for( BoundaryLocation const location : CC::NBS() ) {
      ResetJumpConservatives( location );
   }
}

/**
 * @brief Gives a reference to the corresponding buffer.
 * @param field_type The material field type of the buffer.
 * @param field_index The index of the field asked for.
 * @param conservative_type If a conservative field is wanted, the conservative buffer type. Defaults to ConservativeBufferType::RightHandSide.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetFieldBuffer( MaterialFieldType const field_type, unsigned int const field_index, ConservativeBufferType const conservative_type ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   switch( field_type ) {
      case MaterialFieldType::Conservatives: {
         return GetConservativeBuffer( conservative_type )[field_index];
      }
      case MaterialFieldType::Parameters: {
         return parameters_[field_index];
      }
      default: {// MaterialFieldType::PrimeStates:
         return prime_states_[field_index];
      }
   }
}

/**
 * @brief Const overload.
 */
auto Block::GetFieldBuffer( MaterialFieldType const field_type, unsigned int const field_index, ConservativeBufferType const conservative_type ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   switch( field_type ) {
      case MaterialFieldType::Conservatives: {
         return GetConservativeBuffer( conservative_type )[field_index];
      }
      case MaterialFieldType::Parameters: {
         return parameters_[field_index];
      }
      default: {// MaterialFieldType::PrimeStates:
         return prime_states_[field_index];
      }
   }
}

/**
 * @brief Gives a reference to the corresponding Average buffer.
 * @param equation Decider which buffer is to be returned.
 * @return Reference to the array, that is the requested buffer.
 */
auto Block::GetAverageBuffer( Equation const equation ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return averages_[equation];
}

/**
 * @brief Const overload.
 */
auto Block::GetAverageBuffer( Equation const equation ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return averages_[equation];
}

/**
 * @brief Gives a Reference to the corresponding Right Hand Side buffer.
 * @param equation Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetRightHandSideBuffer( Equation const equation ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return right_hand_sides_[equation];
}

/**
 * @brief Const overload.
 */
auto Block::GetRightHandSideBuffer( Equation const equation ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return right_hand_sides_[equation];
}

/**
 * @brief Gives a Reference to the corresponding initial buffer.
 * @param equation Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetInitialBuffer( Equation const equation ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return initials_[equation];
}

/**
 * @brief Const overload.
 */
auto Block::GetInitialBuffer( Equation const equation ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return initials_[equation];
}

/**
 * @brief Wrapper function that returns a conservative buffer (average, right-hand side or initial). The decision is
 *        made based on the template parameter.
 *        Implementation for the average buffer.
 * @return The average conservative field buffer.
 */
template<>
Conservatives& Block::GetConservativeBuffer<ConservativeBufferType::Average>() {
   return GetAverageBuffer();
}

/**
 * @brief Wrapper function that returns a conservative buffer (average, right-hand side or initial). The decision is
 *        made based on the template parameter.
 *        Implementation for the right-hand side buffer.
 * @return The right-hand side conservative field buffer.
 */
template<>
Conservatives& Block::GetConservativeBuffer<ConservativeBufferType::RightHandSide>() {
   return GetRightHandSideBuffer();
}

/**
 * @brief Wrapper function that returns a conservative buffer (average, right-hand side or initial). The decision is
 *        made based on the template parameter.
 *        Implementation for the initial buffer.
 * @return The initial conservative field buffer.
 */
template<>
Conservatives& Block::GetConservativeBuffer<ConservativeBufferType::Initial>() {
   return GetInitialBuffer();
}

/**
 * @brief Const overload.
 */
template<>
Conservatives const& Block::GetConservativeBuffer<ConservativeBufferType::Average>() const {
   return GetAverageBuffer();
}

/**
 * @brief Const overload.
 */
template<>
Conservatives const& Block::GetConservativeBuffer<ConservativeBufferType::RightHandSide>() const {
   return GetRightHandSideBuffer();
}

/**
 * @brief Const overload.
 */
template<>
Conservatives const& Block::GetConservativeBuffer<ConservativeBufferType::Initial>() const {
   return GetInitialBuffer();
}

/**
 * @brief Gives access to the average buffer.
 * @return Average buffer struct.
 */
Conservatives& Block::GetAverageBuffer() {
   return averages_;
}

/**
 * @brief Const overload.
 */
Conservatives const& Block::GetAverageBuffer() const {
   return averages_;
}

/**
 * @brief Gives access to the right hand side buffer.
 * @return Right-hand side buffer struct.
 */
Conservatives& Block::GetRightHandSideBuffer() {
   return right_hand_sides_;
}

/**
 * @brief Const overload.
 */
Conservatives const& Block::GetRightHandSideBuffer() const {
   return right_hand_sides_;
}

/**
 * @brief Gives access to the initial buffer.
 * @return initial buffer struct.
 */
Conservatives& Block::GetInitialBuffer() {
   return initials_;
}

/**
 * @brief Const overload.
 */
Conservatives const& Block::GetInitialBuffer() const {
   return initials_;
}

/**
 * @brief Gives access to the conservative buffer of given type.
 * @param conservative_type Conservative type of the buffer asked for.
 * @return buffer struct of given type.
 */
Conservatives& Block::GetConservativeBuffer( ConservativeBufferType const conservative_type ) {
   switch( conservative_type ) {
      case ConservativeBufferType::RightHandSide: {
         return right_hand_sides_;
      }
      case ConservativeBufferType::Average: {
         return averages_;
      }
      default: {// case ConservativeBufferType::Initial:
         return initials_;
      }
   }
}

/**
 * @brief Const overload.
 */
Conservatives const& Block::GetConservativeBuffer( ConservativeBufferType const conservative_type ) const {
   switch( conservative_type ) {
      case ConservativeBufferType::RightHandSide: {
         return right_hand_sides_;
      }
      case ConservativeBufferType::Average: {
         return averages_;
      }
      default: {// case ConservativeBufferType::Initial:
         return initials_;
      }
   }
}

/**
 * @brief Gives a Reference to the corresponding Prime-state buffer.
 * @param prime_state_type Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetPrimeStateBuffer( PrimeState const prime_state_type ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return prime_states_[prime_state_type];
}

/**
 * @brief Const overload.
 */
auto Block::GetPrimeStateBuffer( PrimeState const prime_state_type ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return prime_states_[prime_state_type];
}

/**
 * @brief Gives access to the prime state struct.
 * @return Prime state buffer struct.
 */
PrimeStates& Block::GetPrimeStateBuffer() {
   return prime_states_;
}

/**
 * @brief Const overload.
 */
PrimeStates const& Block::GetPrimeStateBuffer() const {
   return prime_states_;
}

/**
 * @brief Gives a Reference to the corresponding Parameter buffer.
 * @param parameter_type Decider which buffer is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetParameterBuffer( Parameter const parameter_type ) -> double ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return parameters_[parameter_type];
}

/**
 * @brief Const overload.
 */
auto Block::GetParameterBuffer( Parameter const parameter_type ) const -> double const ( & )[CC::TCX()][CC::TCY()][CC::TCZ()] {
   return parameters_[parameter_type];
}

/**
 * @brief Gives access to the parameters struct.
 * @return Parameter buffer struct.
 */
Parameters& Block::GetParameterBuffer() {
   return parameters_;
}

/**
 * @brief Const overload.
 */
Parameters const& Block::GetParameterBuffer() const {
   return parameters_;
}

/**
 * @brief Gives a reference to the corresponding jump flux buffer.
 * @param location Decider which face of the block is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetBoundaryJumpFluxes( BoundaryLocation const location ) -> double ( & )[MF::ANOE()][CC::ICY()][CC::ICZ()] {
   return GetBoundaryJump( jump_fluxes_, location );
}

/**
 * @brief Const overload.
 */
auto Block::GetBoundaryJumpFluxes( BoundaryLocation const location ) const -> double const ( & )[MF::ANOE()][CC::ICY()][CC::ICZ()] {
   return GetBoundaryJump( jump_fluxes_, location );
}

/**
 * @brief Gives a reference to the corresponding jump conservative buffer.
 * @param location Decider which face of the block is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetBoundaryJumpConservatives( BoundaryLocation const location ) -> double ( & )[MF::ANOE()][CC::ICY()][CC::ICZ()] {
   return GetBoundaryJump( jump_conservatives_, location );
}

/**
 * @brief Const overload.
 */
auto Block::GetBoundaryJumpConservatives( BoundaryLocation const location ) const -> double const ( & )[MF::ANOE()][CC::ICY()][CC::ICZ()] {
   return GetBoundaryJump( jump_conservatives_, location );
}

/**
 * @brief Gives access to the jump flux buffer.
 * @return The jump flux buffer struct.
 */
SurfaceBuffer& Block::GetBoundaryJumpFluxes() {
   return jump_fluxes_;
}

/**
 * @brief const overload.
 */
SurfaceBuffer const& Block::GetBoundaryJumpFluxes() const {
   return jump_fluxes_;
}

/**
 * @brief Gives access to the jump conservatives buffer.
 * @return The jump conservative buffer struct.
 */
SurfaceBuffer& Block::GetBoundaryJumpConservatives() {
   return jump_conservatives_;
}

/**
 * @brief const overload.
 */
SurfaceBuffer const& Block::GetBoundaryJumpConservatives() const {
   return jump_conservatives_;
}

/**
 * @brief Resets the corresponding jump flux buffer, i.e. set all values in the buffer to zero.
 * @param location Decider which face of the block is to be resetted.
 */
void Block::ResetJumpFluxes( BoundaryLocation const location ) {

   double( &jump )[MF::ANOE()][CC::ICY()][CC::ICZ()] = GetBoundaryJumpFluxes( location );

   for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
      for( unsigned int i = 0; i < CC::ICY(); ++i ) {
         for( unsigned int j = 0; j < CC::ICZ(); ++j ) {
            jump[e][i][j] = 0.0;
         }
      }
   }
}

/**
 * @brief Resets the corresponding jump conservative buffer, i.e. set all values in the buffer to zero.
 * @param location Decider which face of the block is to be resetted.
 */
void Block::ResetJumpConservatives( BoundaryLocation const location ) {

   double( &jump )[MF::ANOE()][CC::ICY()][CC::ICZ()] = GetBoundaryJumpConservatives( location );

   for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
      for( unsigned int i = 0; i < CC::ICY(); ++i ) {
         for( unsigned int j = 0; j < CC::ICZ(); ++j ) {
            jump[e][i][j] = 0.0;
         }
      }
   }
}

/**
 * @brief Gives access to a single conservative array in a SurfaceBuffer struct.
 * @param jump The struct holding the desired array.
 * @param location The location identifer of the desired location.
 * @return Reference to the desired array.
 */
auto GetBoundaryJump( SurfaceBuffer& jump, BoundaryLocation const location ) -> double ( & )[MF::ANOE()][CC::ICY()][CC::ICZ()] {
   switch( location ) {
      case BoundaryLocation::East: {
         return jump.east_;
      }
      case BoundaryLocation::West: {
         return jump.west_;
      }
      case BoundaryLocation::North: {
         return jump.north_;
      }
      case BoundaryLocation::South: {
         return jump.south_;
      }
      case BoundaryLocation::Bottom: {
         return jump.bottom_;
      }
#ifdef PERFORMANCE
      default: {
         return jump.top_;
      }
#else
      case BoundaryLocation::Top: {
         return jump.top_;
      }
      default: {
         throw std::logic_error( "Jump flux Buffer with given Index does not exist" );//suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
      }
#endif
   }
}

/**
 * @brief Const overload.
 */
auto GetBoundaryJump( const SurfaceBuffer& jump, BoundaryLocation const location ) -> double const ( & )[MF::ANOE()][CC::ICY()][CC::ICZ()] {
   switch( location ) {
      case BoundaryLocation::East: {
         return jump.east_;
      }
      case BoundaryLocation::West: {
         return jump.west_;
      }
      case BoundaryLocation::North: {
         return jump.north_;
      }
      case BoundaryLocation::South: {
         return jump.south_;
      }
      case BoundaryLocation::Bottom: {
         return jump.bottom_;
      }
#ifdef PERFORMANCE
      default: {
         return jump.top_;
      }
#else
      case BoundaryLocation::Top: {
         return jump.top_;
      }
      default: {
         throw std::logic_error( "Jump flux Buffer with given Index does not exist" );//suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
      }
#endif
   }
}
