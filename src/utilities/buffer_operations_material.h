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
#ifndef BUFFER_OPERATIONS_MATERIAL_H
#define BUFFER_OPERATIONS_MATERIAL_H

#include "topology/node.h"
#include "utilities/buffer_operations.h"
#include <algorithm>

namespace BufferOperations {

   namespace Material {

      /**
       * @brief This function copies the values in the Conservative buffer from the Source to the Target buffer.
       *        Source and Target ConservativeBufferType are given as templates. It is done for one node.
       * @tparam SourceBuffer The source ConservativeBufferType.
       * @tparam TargetBuffer The target ConservativeBufferType
       * @param node The node for which the buffers are copied.
       */
      template<ConservativeBufferType SourceBuffer, ConservativeBufferType TargetBuffer>
      inline void CopyConservativeBuffersForNode( Node& node ) {
         // Copy operation for all phases
         for( auto& mat_block : node.GetPhases() ) {
            // Get the source and target buffers
            Conservatives const& source_conservatives = mat_block.second.GetConservativeBuffer<SourceBuffer>();
            Conservatives & target_conservatives = mat_block.second.GetConservativeBuffer<TargetBuffer>();
            // Copy the buffer
            BO::CopyFieldBuffer( source_conservatives, target_conservatives );
         }
      }

      /**
       * @brief This function copies the values in the Conservative buffer from the Source to the Target buffer.
       *        Source and Target ConservativeBufferType are given as templates. It is done for a vector of nodes. Calls a copy
       *        function for each node in the vector.
       * @tparam SourceBuffer The source ConservativeBufferType.
       * @tparam TargetBuffer The target ConservativeBufferType
       * @param nodes The nodes for which the buffers are copied.
       */
      template<ConservativeBufferType SourceBuffer, ConservativeBufferType TargetBuffer>
      inline void CopyConservativeBuffersForNodeList( std::vector<std::reference_wrapper<Node>> const& nodes ) {
         for( Node& node : nodes ) {
            CopyConservativeBuffersForNode<SourceBuffer, TargetBuffer>( node );
         }
      }

      /**
       * @brief Swaps the Conservative buffer of the FirstBuffer and SecondBuffer ConservativeBufferType.
       *        This is done for a single node.
       * @tparam FirstBuffer The first ConservativeBufferType.
       * @tparam SecondBuffer The second ConservativeBufferType.
       * @param node The node for which the buffers are swapped.
       */
      template<ConservativeBufferType FirstBuffer, ConservativeBufferType SecondBuffer>
      inline void SwapConservativeBuffersForNode( Node & node ) {
         for( auto& mat_block : node.GetPhases() ) {
            Conservatives & first_conservatives  = mat_block.second.GetConservativeBuffer<FirstBuffer>();
            Conservatives & second_conservatives = mat_block.second.GetConservativeBuffer<SecondBuffer>();
            BO::SwapFieldBuffer( first_conservatives, second_conservatives );
         }
      }

      /**
       * @brief Swaps the Conservative buffer of the FirstBuffer and SecondBuffer ConservativeBufferType.
       *        This is done for a node vector. Calls the a swap function on each node in the vector.
       * @tparam FirstBuffer The first ConservativeBufferType.
       * @tparam SecondBuffer The second ConservativeBufferType.
       * @param node The nodes for which the buffers are swapped.
       */
      template<ConservativeBufferType FirstBuffer, ConservativeBufferType SecondBuffer>
      inline void SwapConservativeBuffersForNodeList( std::vector<std::reference_wrapper<Node>> const& nodes ) {
         for( Node & node : nodes ) {
            SwapConservativeBuffersForNode<FirstBuffer, SecondBuffer>( node );
         }
      }

   } // namespace Material

} // namespace BufferOperations

#endif // BUFFER_OPERATIONS_MATERIAL_H