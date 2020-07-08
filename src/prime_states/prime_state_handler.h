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
#ifndef PRIME_STATE_HANDLER_H
#define PRIME_STATE_HANDLER_H

#include "block_definitions/block.h"
#include "block_definitions/field_material_definitions.h"
#include "materials/material_manager.h"

/**
 * @brief This class handles the conversion between conservative and prime state values for full buffers or single cells.
 * @tparam DerivedPrimeStateHandler CRTP template parameter.
 */
template<typename DerivedPrimeStateHandler>
class PrimeStateHandler {

   friend DerivedPrimeStateHandler;

   MaterialManager const& material_manager_;

   explicit PrimeStateHandler( MaterialManager const& material_manager ) : material_manager_( material_manager ) {}

public:
   ~PrimeStateHandler()                          = default;
   PrimeStateHandler( PrimeStateHandler const& ) = delete;
   PrimeStateHandler& operator=( PrimeStateHandler const& ) = delete;
   PrimeStateHandler( PrimeStateHandler&& )                 = delete;
   PrimeStateHandler& operator=( PrimeStateHandler&& ) = delete;

   /**
    * @brief Converts prime state to conservative values for every cell in the given buffers.
    * @param material Material identifier of the material under consideration.
    * @param prime_states The input PrimeStates buffer.
    * @param conservatives The output Conservatives buffer.
    */
   void ConvertPrimeStatesToConservatives( MaterialName const& material, PrimeStates const& prime_states, Conservatives& conservatives ) const {
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               ConvertPrimeStatesToConservatives( material, prime_states, conservatives, i, j, k );
            }
         }
      }
   }

   /**
    * @brief Converts prime state to conservative values for a specific cell in the given buffers.
    * @param material Material identifier of the material under consideration.
    * @param prime_states The input PrimeStates buffer.
    * @param conservatives The output Conservatives buffer.
    * @param i,j,k Indices of the cell.
    */
   inline void ConvertPrimeStatesToConservatives( MaterialName const& material, PrimeStates const& prime_states, Conservatives& conservatives, unsigned int const i, unsigned int const j, unsigned int const k ) const {
      auto&& conservatives_cell = conservatives.GetCellView( i, j, k );
      ConvertPrimeStatesToConservatives( material, prime_states.GetCellView( i, j, k ), conservatives_cell );
   }

   /**
    * @brief Converts prime state to conservative values.
    * @param material Material identifier of the material under consideration.
    * @param prime_states_container Container of the input prime states.
    * @param conservatives_container Container of the output conservatives.
    * @tparam PrimeStatesContainerType Type of the conservative container, has to provide an [index]-operator.
    * @tparam ConservativesContainerType Type of the prime state container, has to provide an [index]-operator.
    * @note In the containers, prime state and conservative values are assumed to be at the indices given by the PrimeState and Equation enumerations, respectively.
    */
   template<typename PrimeStatesContainerType, typename ConservativesContainerType>
   inline void ConvertPrimeStatesToConservatives( MaterialName const& material, PrimeStatesContainerType const& prime_states_container, ConservativesContainerType& conservatives_container ) const {
      static_cast<DerivedPrimeStateHandler const&>( *this ).ConvertPrimeStatesToConservativesImplementation( material, prime_states_container, conservatives_container );
   }

   /**
    * @brief Converts conservative to prime state values for every cell in the given buffers.
    * @param material Material identifier of the material under consideration.
    * @param conservatives The input Conservatives buffer.
    * @param prime_states The output PrimeStates buffer.
    */
   void ConvertConservativesToPrimeStates( MaterialName const& material, Conservatives const& conservatives, PrimeStates& prime_states ) const {
      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               ConvertConservativesToPrimeStates( material, conservatives, prime_states, i, j, k );
            }
         }
      }
   }

   /**
    * @brief Converts conservative to prime state values for a specific cell in the given buffers.
    * @param material Material identifier of the material under consideration.
    * @param conservatives The input Conservatives buffer.
    * @param prime_states The output PrimeStates buffer.
    * @param i,j,k Indices of the cell.
    */
   inline void ConvertConservativesToPrimeStates( MaterialName const& material, Conservatives const& conservatives, PrimeStates& prime_states, unsigned int const i, unsigned int const j, unsigned int const k ) const {
      auto&& prime_states_cell = prime_states.GetCellView( i, j, k );
      ConvertConservativesToPrimeStates( material, conservatives.GetCellView( i, j, k ), prime_states_cell );
   }

   /**
    * @brief Converts conservative to prime state values.
    * @param material Material identifier of the material under consideration.
    * @param conservatives_container Container of the input conservatives.
    * @param prime_states_container Container of the output prime states.
    * @tparam ConservativesContainerType Type of the prime state container, has to provide an [index]-operator.
    * @tparam PrimeStatesContainerType Type of the conservative container, has to provide an [index]-operator.
    * @note In the containers, prime state and conservative values are assumed to be at the indices given by the PrimeState and Equation enumerations, respectively.
    */
   template<typename ConservativesContainerType, typename PrimeStatesContainerType>
   inline void ConvertConservativesToPrimeStates( MaterialName const& material, ConservativesContainerType const& conservatives_container, PrimeStatesContainerType& prime_states_container ) const {
      static_cast<DerivedPrimeStateHandler const&>( *this ).ConvertConservativesToPrimeStatesImplementation( material, conservatives_container, prime_states_container );
   }
};

#endif// PRIME_STATE_HANDLER_H
