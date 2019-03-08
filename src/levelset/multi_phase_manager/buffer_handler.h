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
#ifndef BUFFER_HANDLER_H
#define BUFFER_HANDLER_H


#include "materials/material_manager.h"
#include "topology/node.h"

/**
 * @brief The BufferHandler class provides functionality to translate values from the different buffers (e.g. calculate prime state from conservatives for
 * only real fluid cells).
 * @tparam Typename as template parameter due to CRTP.
 */
template<typename DerivedBufferHandler>
class BufferHandler {

   friend DerivedBufferHandler;

   MaterialManager const& material_manager_;

   /**
    * @brief Constructor for the buffer handler for level-set simulations.
    * @param material_manager Instance of a material manager, which already has been initialized according to the user input.
    */
   explicit BufferHandler( MaterialManager const& material_manager) : material_manager_(material_manager) { }

public:
   BufferHandler() = delete;
   ~BufferHandler() = default;
   BufferHandler( BufferHandler const& ) = delete;
   BufferHandler& operator=( BufferHandler const& )= delete;
   BufferHandler( BufferHandler&& ) = delete;
   BufferHandler& operator=( BufferHandler&& ) = delete;

   /**
    * @brief Transform given volume averaged conservatives to conservatives. This is done by a multiplication with the volume fraction.
    * @param node The node for which conservatives are calculated.
    */
   void TransformToConservatives( Node& node ) const {
      static_cast<DerivedBufferHandler const&>( *this ).TransformToConservativesImplementation( node );
   }

   /**
    * @brief Transform given conservatives to volume averaged conservatives. This is done by a division with the volume fraction.
    * @param node The node for which volume averaged conservatives are calculated.
    */
   void TransformToVolumeAveragedConservatives( Node& node ) const {
      static_cast<DerivedBufferHandler const&>( *this ).TransformToVolumeAveragedConservativesImplementation( node );
   }

   /**
    * @brief During the scale-separation procedure small flow structures at the interface get dissolved. Thus, in the material which is not dissolved,
    * real-material cells can be generated. Those cells have to be filled with prime-state values from the last RK stage.
    * @param node The node, for which the conservatives have to be corrected.
    */
   void AdaptConservativesToWellResolvedDistanceFunction( Node& node ) const {
      static_cast<DerivedBufferHandler const&>( *this ).AdaptConservativesToWellResolvedDistanceFunctionImplementation( node );
   }

   /**
    * @brief We integrate conservatives in time. After time integration it is necessary to calculate and store the prime states
    * for the integrated conservatives. This is done in this function.
    * @param node The node for which we calculate the prime states.
    */
   void CalculatePrimesFromIntegratedConservatives( Node& node ) const {
      static_cast<DerivedBufferHandler const&>( *this ).CalculatePrimesFromIntegratedConservativesImplementation( node) ;
   }

   /**
    * @brief Populates the cells of the conservative_rhs in which we extendwith correct values, based on the information we have in the
    * prime state buffer.
    * @param node The node for which conservatives are calculated.
    */
   void CalculateConservativesFromExtendedPrimes( Node& node ) const {
      static_cast<DerivedBufferHandler const&>( *this ).CalculateConservativesFromExtendedPrimesImplementation( node );
   }
};


#endif //BUFFER_HANDLER_H
