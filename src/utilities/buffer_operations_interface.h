//===----------------- buffer_operations_interface.h ----------------------===//
//
//                                 ALPACA
//
// Part of ALPACA, under the GNU General Public License as published by
// the Free Software Foundation version 3.
// SPDX-License-Identifier: GPL-3.0-only
//
// If using this code in an academic setting, please cite the following:
// @article{hoppe2022parallel,
//  title={A parallel modular computing environment for three-dimensional
//  multiresolution simulations of compressible flows},
//  author={Hoppe, Nils and Adami, Stefan and Adams, Nikolaus A},
//  journal={Computer Methods in Applied Mechanics and Engineering},
//  volume={391},
//  pages={114486},
//  year={2022},
//  publisher={Elsevier}
// }
//
//===----------------------------------------------------------------------===//
#ifndef BUFFER_OPERATIONS_INTERFACE_H
#define BUFFER_OPERATIONS_INTERFACE_H

#include "topology/node.h"
#include "utilities/buffer_operations.h"
#include <algorithm>

namespace BufferOperations {

namespace Interface {

/**
 * @brief This function copies the values in the InterfaceDescription buffer
 * from the Source to the Target buffer. Source and Target
 * InterfaceDescriptionBufferType are given as templates. It is done for one
 * node.
 * @tparam SourceBuffer The source InterfaceDescriptionBufferType.
 * @tparam TargetBuffer The target InterfaceDescriptionBufferType.
 * @tparam The InterfaceDescription type that should be copied (Default:
 * Levelset).
 * @param node The node for which the buffers are copied.
 */
template <InterfaceDescriptionBufferType SourceBuffer,
          InterfaceDescriptionBufferType TargetBuffer,
          InterfaceDescription Type = InterfaceDescription::Levelset>
inline void CopyInterfaceDescriptionBufferForNode(Node &node) {
  // Get the source and target buffers
  InterfaceBlock &interface_block = node.GetInterfaceBlock();
  double const(&source_description)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetInterfaceDescriptionBuffer<SourceBuffer>()[Type];
  double(&target_description)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      interface_block.GetInterfaceDescriptionBuffer<TargetBuffer>()[Type];
  // Copy the values
  BO::CopySingleBuffer(source_description, target_description);
}

/**
 * @brief This function copies the values in the InterfaceDescription buffer
 * from the Source to the Target buffer. Source and Target
 * InterfaceDescriptionBufferType are given as templates. It is done for a
 * vector of nodes. Calls a copy function for each node in the vector.
 * @tparam SourceBuffer The source InterfaceDescriptionBufferType.
 * @tparam TargetBuffer The target InterfaceDescriptionBufferType.
 * @tparam The InterfaceDescription type that should be copied (Default:
 * Levelset).
 * @param nodes The nodes for which the buffers are copied.
 */
template <InterfaceDescriptionBufferType SourceBuffer,
          InterfaceDescriptionBufferType TargetBuffer,
          InterfaceDescription Type = InterfaceDescription::Levelset>
inline void CopyInterfaceDescriptionBufferForNodeList(
    std::vector<std::reference_wrapper<Node>> const &nodes) {
  for (Node &node : nodes) {
    CopyInterfaceDescriptionBufferForNode<SourceBuffer, TargetBuffer, Type>(
        node);
  }
}

/**
 * @brief Swaps the InterfaceDescription buffer of the FirstBuffer and
 * SecondBuffer InterfaceDescriptionBufferType. This is done for a single node.
 * @tparam FirstBuffer The first InterfaceDescriptionBufferType.
 * @tparam SecondBuffer The second InterfaceDescriptionBufferType.
 * @tparam The InterfaceDescription type that should be copied (Default:
 * Levelset).
 * @param node The node for which the buffers are swapped.
 */
template <InterfaceDescriptionBufferType FirstBuffer,
          InterfaceDescriptionBufferType SecondBuffer,
          InterfaceDescription Type = InterfaceDescription::Levelset>
inline void SwapInterfaceDescriptionBufferForNode(Node &node) {
  double(&first_description_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock()
          .GetInterfaceDescriptionBuffer<FirstBuffer>()[Type];
  double(&second_description_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
      node.GetInterfaceBlock()
          .GetInterfaceDescriptionBuffer<SecondBuffer>()[Type];
  BO::SwapSingleBuffer(first_description_buffer, second_description_buffer);
}

/**
 * @brief Swaps the InterfaceDescription buffer of the FirstBuffer and
 * SecondBuffer InterfaceDescriptionBufferType. This is done for a node vector.
 * Calls the a swap function on each node in the vector.
 * @tparam FirstBuffer The first InterfaceDescriptionBufferType.
 * @tparam SecondBuffer The second InterfaceDescriptionBufferType.
 * @tparam The InterfaceDescription type that should be copied (Default:
 * Levelset).
 * @param node The nodes for which the buffers are swapped.
 */
template <InterfaceDescriptionBufferType FirstBuffer,
          InterfaceDescriptionBufferType SecondBuffer,
          InterfaceDescription Type = InterfaceDescription::Levelset>
inline void SwapInterfaceDescriptionBufferForNodeList(
    std::vector<std::reference_wrapper<Node>> const &nodes) {
  for (Node &node : nodes) {
    SwapInterfaceDescriptionBufferForNode<FirstBuffer, SecondBuffer, Type>(
        node);
  }
}

} // namespace Interface

} // namespace BufferOperations

#endif // BUFFER_OPERATIONS_INTERFACE_H
