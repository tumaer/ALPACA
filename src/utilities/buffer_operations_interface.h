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
#ifndef BUFFER_OPERATIONS_INTERFACE_H
#define BUFFER_OPERATIONS_INTERFACE_H

#include "topology/node.h"
#include "utilities/buffer_operations.h"
#include <algorithm>

namespace BufferOperations {

   namespace Interface {

      /**
       * @brief This function copies the values in the InterfaceDescription buffer from the Source to the Target buffer.
       *        Source and Target InterfaceDescriptionBufferType are given as templates. It is done for one node.
       * @tparam SourceBuffer The source InterfaceDescriptionBufferType.
       * @tparam TargetBuffer The target InterfaceDescriptionBufferType.
       * @tparam The InterfaceDescription type that should be copied (Default: Levelset).
       * @param node The node for which the buffers are copied.
       */
      template<InterfaceDescriptionBufferType SourceBuffer, InterfaceDescriptionBufferType TargetBuffer, InterfaceDescription Type = InterfaceDescription::Levelset>
      inline void CopyInterfaceDescriptionBufferForNode( Node& node ) {
         // Get the source and target buffers
         InterfaceBlock& interface_block = node.GetInterfaceBlock();
         double const (&source_description)[CC::TCX()][CC::TCY()][CC::TCZ()] = interface_block.GetInterfaceDescriptionBuffer<SourceBuffer>()[Type];
         double (&target_description)[CC::TCX()][CC::TCY()][CC::TCZ()] = interface_block.GetInterfaceDescriptionBuffer<TargetBuffer>()[Type];
         // Copy the values
         BO::CopySingleBuffer( source_description, target_description );
      }

      /**
       * @brief This function copies the values in the InterfaceDescription buffer from the Source to the Target buffer.
       *        Source and Target InterfaceDescriptionBufferType are given as templates. It is done for a vector of nodes. Calls a copy
       *        function for each node in the vector.
       * @tparam SourceBuffer The source InterfaceDescriptionBufferType.
       * @tparam TargetBuffer The target InterfaceDescriptionBufferType.
       * @tparam The InterfaceDescription type that should be copied (Default: Levelset).
       * @param nodes The nodes for which the buffers are copied.
       */
      template<InterfaceDescriptionBufferType SourceBuffer, InterfaceDescriptionBufferType TargetBuffer, InterfaceDescription Type = InterfaceDescription::Levelset>
      inline void CopyInterfaceDescriptionBufferForNodeList( std::vector<std::reference_wrapper<Node>> const& nodes ) {
         for( Node& node : nodes ) {
            CopyInterfaceDescriptionBufferForNode<SourceBuffer, TargetBuffer, Type>( node );
         }
      }

      /**
       * @brief Swaps the InterfaceDescription buffer of the FirstBuffer and SecondBuffer InterfaceDescriptionBufferType.
       *        This is done for a single node.
       * @tparam FirstBuffer The first InterfaceDescriptionBufferType.
       * @tparam SecondBuffer The second InterfaceDescriptionBufferType.
       * @tparam The InterfaceDescription type that should be copied (Default: Levelset).
       * @param node The node for which the buffers are swapped.
       */
      template<InterfaceDescriptionBufferType FirstBuffer, InterfaceDescriptionBufferType SecondBuffer, InterfaceDescription Type = InterfaceDescription::Levelset>
      inline void SwapInterfaceDescriptionBufferForNode( Node & node ) {
         double (&first_description_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetInterfaceDescriptionBuffer<FirstBuffer>()[Type];
         double (&second_description_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetInterfaceDescriptionBuffer<SecondBuffer>()[Type];
         BO::SwapSingleBuffer( first_description_buffer, second_description_buffer );
      }

      /**
       * @brief Swaps the InterfaceDescription buffer of the FirstBuffer and SecondBuffer InterfaceDescriptionBufferType.
       *        This is done for a node vector. Calls the a swap function on each node in the vector.
       * @tparam FirstBuffer The first InterfaceDescriptionBufferType.
       * @tparam SecondBuffer The second InterfaceDescriptionBufferType.
       * @tparam The InterfaceDescription type that should be copied (Default: Levelset).
       * @param node The nodes for which the buffers are swapped.
       */
      template<InterfaceDescriptionBufferType FirstBuffer, InterfaceDescriptionBufferType SecondBuffer, InterfaceDescription Type = InterfaceDescription::Levelset>
      inline void SwapInterfaceDescriptionBufferForNodeList( std::vector<std::reference_wrapper<Node>> const& nodes ) {
         for( Node & node : nodes ) {
            SwapInterfaceDescriptionBufferForNode<FirstBuffer, SecondBuffer, Type>( node );
         }
      }

   } // namespace Interface

} // namespace BufferOperations

#endif // BUFFER_OPERATIONS_INTERFACE_H