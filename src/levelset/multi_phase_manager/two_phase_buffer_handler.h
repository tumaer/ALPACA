//===-------------------- two_phase_buffer_handler.h ----------------------===//
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
#ifndef TWO_PHASE_BUFFER_HANDLER_H
#define TWO_PHASE_BUFFER_HANDLER_H

#include "buffer_handler.h"
#include "prime_states/prime_state_handler.h"

/**
 * @brief The TwoPhaseBufferHandler class is a specification of the
 * BufferHandler for two phase simulations.
 */
class TwoPhaseBufferHandler : public BufferHandler<TwoPhaseBufferHandler> {

  friend BufferHandler;

  PrimeStateHandler const prime_state_handler_;

private:
  void TransformToConservativesImplementation(Node &node) const;
  void TransformToVolumeAveragedConservativesImplementation(Node &node) const;
  void AdaptConservativesToWellResolvedDistanceFunctionImplementation(
      Node &node) const;
  void
  CalculatePrimesFromIntegratedConservativesImplementation(Node &node) const;
  void CalculateConservativesFromExtendedPrimesImplementation(Node &node) const;

public:
  TwoPhaseBufferHandler() = delete;
  explicit TwoPhaseBufferHandler(MaterialManager const &material_manager);
  ~TwoPhaseBufferHandler() = default;
  TwoPhaseBufferHandler(TwoPhaseBufferHandler const &) = delete;
  TwoPhaseBufferHandler &operator=(TwoPhaseBufferHandler const &) = delete;
  TwoPhaseBufferHandler(TwoPhaseBufferHandler &&) = delete;
  TwoPhaseBufferHandler &operator=(TwoPhaseBufferHandler &&) = delete;
};

#endif // TWO_PHASE_BUFFER_HANDLER_H
