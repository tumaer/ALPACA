//===----------------- buffer_operations_material.h -----------------------===//
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
#ifndef BUFFER_OPERATIONS_MATERIAL_H
#define BUFFER_OPERATIONS_MATERIAL_H

#include "topology/node.h"
#include "utilities/buffer_operations.h"
#include <algorithm>

namespace BufferOperations {

namespace Material {

/**
 * @brief This function copies the values in the Conservative buffer from the
 * Source to the Target buffer. Source and Target ConservativeBufferType are
 * given as templates. It is done for one node.
 * @tparam SourceBuffer The source ConservativeBufferType.
 * @tparam TargetBuffer The target ConservativeBufferType.
 * @param node The node for which the buffers are copied.
 */
template <ConservativeBufferType SourceBuffer,
          ConservativeBufferType TargetBuffer>
inline void CopyConservativeBuffersForNode(Node &node) {
  // Copy operation for all phases
  for (auto &mat_block : node.GetPhases()) {
    // Get the source and target buffers
    Conservatives const &source_conservatives =
        mat_block.second.GetConservativeBuffer<SourceBuffer>();
    Conservatives &target_conservatives =
        mat_block.second.GetConservativeBuffer<TargetBuffer>();
    // Copy the buffer
    BO::CopyFieldBuffer(source_conservatives, target_conservatives);
  }
}

/**
 * @brief This function copies the values in the Conservative buffer from the
 * Source to the Target buffer. Source and Target ConservativeBufferType are
 * given as templates. It is done for a vector of nodes. Calls a copy function
 * for each node in the vector.
 * @tparam SourceBuffer The source ConservativeBufferType.
 * @tparam TargetBuffer The target ConservativeBufferType.
 * @param nodes The nodes for which the buffers are copied.
 */
template <ConservativeBufferType SourceBuffer,
          ConservativeBufferType TargetBuffer>
inline void CopyConservativeBuffersForNodeList(
    std::vector<std::reference_wrapper<Node>> const &nodes) {
  for (Node &node : nodes) {
    CopyConservativeBuffersForNode<SourceBuffer, TargetBuffer>(node);
  }
}

/**
 * @brief Swaps the Conservative buffer of the FirstBuffer and SecondBuffer
 * ConservativeBufferType. This is done for a single node.
 * @tparam FirstBuffer The first ConservativeBufferType.
 * @tparam SecondBuffer The second ConservativeBufferType.
 * @param node The node for which the buffers are swapped.
 */
template <ConservativeBufferType FirstBuffer,
          ConservativeBufferType SecondBuffer>
inline void SwapConservativeBuffersForNode(Node &node) {
  for (auto &mat_block : node.GetPhases()) {
    Conservatives &first_conservatives =
        mat_block.second.GetConservativeBuffer<FirstBuffer>();
    Conservatives &second_conservatives =
        mat_block.second.GetConservativeBuffer<SecondBuffer>();
    BO::SwapFieldBuffer(first_conservatives, second_conservatives);
  }
}

/**
 * @brief Swaps the Conservative buffer of the FirstBuffer and SecondBuffer
 * ConservativeBufferType. This is done for a node vector. Calls the a swap
 * function on each node in the vector.
 * @tparam FirstBuffer The first ConservativeBufferType.
 * @tparam SecondBuffer The second ConservativeBufferType.
 * @param node The nodes for which the buffers are swapped.
 */
template <ConservativeBufferType FirstBuffer,
          ConservativeBufferType SecondBuffer>
inline void SwapConservativeBuffersForNodeList(
    std::vector<std::reference_wrapper<Node>> const &nodes) {
  for (Node &node : nodes) {
    SwapConservativeBuffersForNode<FirstBuffer, SecondBuffer>(node);
  }
}

} // namespace Material

} // namespace BufferOperations

#endif // BUFFER_OPERATIONS_MATERIAL_H
