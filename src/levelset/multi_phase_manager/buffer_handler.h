//===------------------------ buffer_handler.h ----------------------------===//
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
#ifndef BUFFER_HANDLER_H
#define BUFFER_HANDLER_H

#include "materials/material_manager.h"
#include "topology/node.h"

/**
 * @brief The BufferHandler class provides functionality to translate values
 * from the different buffers (e.g. calculate prime state from conservatives for
 *        only real fluid cells).
 * @tparam Typename as template parameter due to CRTP.
 */
template <typename DerivedBufferHandler> class BufferHandler {

  friend DerivedBufferHandler;

  MaterialManager const &material_manager_;

  /**
   * @brief Constructor for the buffer handler for level-set simulations.
   * @param material_manager Instance of a material manager, which already has
   * been initialized according to the user input.
   */
  explicit BufferHandler(MaterialManager const &material_manager)
      : material_manager_(material_manager) {}

public:
  BufferHandler() = delete;
  ~BufferHandler() = default;
  BufferHandler(BufferHandler const &) = delete;
  BufferHandler &operator=(BufferHandler const &) = delete;
  BufferHandler(BufferHandler &&) = delete;
  BufferHandler &operator=(BufferHandler &&) = delete;

  /**
   * @brief Transform given volume averaged conservatives to conservatives. This
   * is done by a multiplication with the volume fraction.
   * @param node The node for which conservatives are calculated.
   */
  void TransformToConservatives(Node &node) const {
    static_cast<DerivedBufferHandler const &>(*this)
        .TransformToConservativesImplementation(node);
  }

  /**
   * @brief Transform given conservatives to volume averaged conservatives. This
   * is done by a division with the volume fraction.
   * @param node The node for which volume averaged conservatives are
   * calculated.
   */
  void TransformToVolumeAveragedConservatives(Node &node) const {
    static_cast<DerivedBufferHandler const &>(*this)
        .TransformToVolumeAveragedConservativesImplementation(node);
  }

  /**
   * @brief During the scale-separation procedure small flow structures at the
   * interface get dissolved. Thus, in the material which is not dissolved,
   *        real-material cells can be generated. Those cells have to be filled
   * with prime-state values from the last RK stage.
   * @param node The node, for which the conservatives have to be corrected.
   */
  void AdaptConservativesToWellResolvedDistanceFunction(Node &node) const {
    static_cast<DerivedBufferHandler const &>(*this)
        .AdaptConservativesToWellResolvedDistanceFunctionImplementation(node);
  }

  /**
   * @brief We integrate conservatives in time. After time integration it is
   * necessary to calculate and store the prime states for the integrated
   * conservatives. This is done in this function.
   * @param node The node for which we calculate the prime states.
   */
  void CalculatePrimesFromIntegratedConservatives(Node &node) const {
    static_cast<DerivedBufferHandler const &>(*this)
        .CalculatePrimesFromIntegratedConservativesImplementation(node);
  }

  /**
   * @brief Populates the cells of the conservative_rhs in which we extendwith
   * correct values, based on the information we have in the prime state buffer.
   * @param node The node for which conservatives are calculated.
   */
  void CalculateConservativesFromExtendedPrimes(Node &node) const {
    static_cast<DerivedBufferHandler const &>(*this)
        .CalculateConservativesFromExtendedPrimesImplementation(node);
  }
};

#endif // BUFFER_HANDLER_H
