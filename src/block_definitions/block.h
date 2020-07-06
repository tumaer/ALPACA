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
#ifndef BLOCK_H
#define BLOCK_H

#include "user_specifications/compile_time_constants.h"
#include "boundary_condition/boundary_specifications.h"
#include "block_definitions/field_buffer.h"
#include "block_definitions/field_material_definitions.h"

/**
 * @brief Gives a buffer for the values on the six (in 3D) surfaces of the block.
 */
struct SurfaceBuffer{
   double east_[MF::ANOE()][CC::ICY()][CC::ICZ()];
   double west_[MF::ANOE()][CC::ICY()][CC::ICZ()];
   double north_[MF::ANOE()][CC::ICY()][CC::ICZ()];
   double south_[MF::ANOE()][CC::ICY()][CC::ICZ()];
   double top_[MF::ANOE()][CC::ICY()][CC::ICZ()];
   double bottom_[MF::ANOE()][CC::ICY()][CC::ICZ()];
};
// Check Memory Layout at compile time for safe MPI sending (Ensures Compiler did not pad the struct)
static_assert( sizeof(SurfaceBuffer) == 6 * MF::ANOE() * CC::ICY() * CC::ICZ() * sizeof(double), "Surface Struct is not contiguous in Memory" );

/**
 * @brief The Block class holds the data on which the simulation is running. They do NOT manipulate the data themselves, but provide data access
 *        to the different field buffers. One Block only contains material data for one material. >>A block is always single-phase<<.
 */
class Block {
   // buffers for the conservatives (different buffer types required for the integration)
   Conservatives averages_;
   Conservatives right_hand_sides_;
   Conservatives initials_;

   // buffers for the primestates (e.g. temperature, pressure, velocity)
   PrimeStates prime_states_;

   //buffer to store location-dependent material parameter (e.g., viscosity, conductivity)
   Parameters parameters_;

   //buffer to save fluxes at internal jump boundaries
   SurfaceBuffer jump_fluxes_;

   //buffer to store conservative fluxes at internal jump boundaries
   SurfaceBuffer jump_conservatives_;

public:
   explicit Block();
   ~Block() = default;
   Block( Block const& ) = delete;
   Block& operator=( Block const& ) = delete;
   Block( Block&& ) = delete;
   Block& operator=( Block&& ) = delete;

   // Returning general field buffer
   auto GetFieldBuffer( MaterialFieldType const field_type, unsigned int const field_index, ConservativeBufferType const conservative_type = ConservativeBufferType::RightHandSide ) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetFieldBuffer( MaterialFieldType const field_type, unsigned int const field_index, ConservativeBufferType const conservative_type = ConservativeBufferType::RightHandSide ) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   // Returning conservative buffers
   auto GetAverageBuffer( Equation const equation ) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetAverageBuffer( Equation const equation ) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetRightHandSideBuffer( Equation const equation ) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetRightHandSideBuffer( Equation const equation ) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   auto GetInitialBuffer( Equation const equation ) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetInitialBuffer( Equation const equation ) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   template<ConservativeBufferType C>
   Conservatives& GetConservativeBuffer();

   template<ConservativeBufferType C>
   Conservatives const& GetConservativeBuffer() const;

   Conservatives& GetAverageBuffer();
   Conservatives const& GetAverageBuffer() const;
   Conservatives& GetRightHandSideBuffer();
   Conservatives const& GetRightHandSideBuffer() const;
   Conservatives& GetInitialBuffer();
   Conservatives const& GetInitialBuffer() const;
   Conservatives& GetConservativeBuffer( ConservativeBufferType const conservative_type );
   Conservatives const& GetConservativeBuffer( ConservativeBufferType const conservative_type ) const;

   // Returning primestate buffers
   auto GetPrimeStateBuffer( PrimeState const prime_state_type ) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetPrimeStateBuffer( PrimeState const prime_state_type ) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   PrimeStates& GetPrimeStateBuffer();
   PrimeStates const& GetPrimeStateBuffer() const;

   // returning parameter buffers
   auto GetParameterBuffer( Parameter const parameter_type ) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
   auto GetParameterBuffer( Parameter const parameter_type ) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

   Parameters& GetParameterBuffer();
   Parameters const& GetParameterBuffer() const;

   // Returning surface buffer
   auto GetBoundaryJumpFluxes( BoundaryLocation const location ) -> double (&)[MF::ANOE()][CC::ICY()][CC::ICZ()];
   auto GetBoundaryJumpFluxes( BoundaryLocation const location ) const -> double const (&)[MF::ANOE()][CC::ICY()][CC::ICZ()];
   auto GetBoundaryJumpConservatives( BoundaryLocation const location ) -> double (&)[MF::ANOE()][CC::ICY()][CC::ICZ()];
   auto GetBoundaryJumpConservatives( BoundaryLocation const location ) const -> double const (&)[MF::ANOE()][CC::ICY()][CC::ICZ()];

   SurfaceBuffer& GetBoundaryJumpFluxes();
   SurfaceBuffer const& GetBoundaryJumpFluxes() const;
   SurfaceBuffer& GetBoundaryJumpConservatives();
   SurfaceBuffer const& GetBoundaryJumpConservatives() const;

   void ResetJumpFluxes( BoundaryLocation const location );
   void ResetJumpConservatives( BoundaryLocation const location );
};

auto GetBoundaryJump( SurfaceBuffer& jump, BoundaryLocation const location ) -> double (&)[MF::ANOE()][CC::ICY()][CC::ICZ()];
auto GetBoundaryJump( SurfaceBuffer const& jump, BoundaryLocation const location ) -> double const (&)[MF::ANOE()][CC::ICY()][CC::ICZ()];

#endif // BLOCK_H