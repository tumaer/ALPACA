//===------------------------ interface_block.h ---------------------------===//
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
#ifndef INTERFACE_BLOCK_H
#define INTERFACE_BLOCK_H

#include "block_definitions/field_buffer.h"
#include "block_definitions/field_interface_definitions.h"
#include "interface_block_buffer_definitions.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief The InterfaceBlock class holds the interface data, such as the
 * interface describing variables (levelset, volume fraction) and other
 * interfaces variables that are required for the computations (e.g., interface
 * velocity) Does NOT manipulate the data itself, but provides access to the
 * data.
 */
class InterfaceBlock {

  // buffers for the interface description (different buffer types required for
  // the integration)
  InterfaceDescriptions base_;
  InterfaceDescriptions right_hand_side_;
  InterfaceDescriptions reinitialized_;
  InterfaceDescriptions initial_;
  InterfaceDescriptions integrated_;

  // buffers for the interface states (e.g. interface velocity,
  // negative/positive pressure)
  InterfaceStates states_;

  // buffer to store field-dependent interface parameter (e.g., surface tension
  // coefficient)
  InterfaceParameters parameters_;

public:
  InterfaceBlock() = delete;
  explicit InterfaceBlock(double const levelset_initial);
  explicit InterfaceBlock(
      double const (&levelset_initial)[CC::TCX()][CC::TCY()][CC::TCZ()]);
  ~InterfaceBlock() = default;
  InterfaceBlock(InterfaceBlock const &) = delete;
  InterfaceBlock &operator=(InterfaceBlock const &) = delete;
  InterfaceBlock(InterfaceBlock &&) = delete;
  InterfaceBlock &operator=(InterfaceBlock &&) = delete;

  // Returning general field buffer
  auto GetFieldBuffer(InterfaceFieldType const field_type,
                      unsigned int const field_index,
                      InterfaceDescriptionBufferType const buffer_type =
                          InterfaceDescriptionBufferType::Base)
      -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto GetFieldBuffer(InterfaceFieldType const field_type,
                      unsigned int const field_index,
                      InterfaceDescriptionBufferType const buffer_type =
                          InterfaceDescriptionBufferType::Base) const
      -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

  // Returning interface description buffers
  auto GetBaseBuffer(InterfaceDescription const interface_description)
      -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto GetBaseBuffer(InterfaceDescription const interface_description) const
      -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

  auto GetRightHandSideBuffer(InterfaceDescription const interface_description)
      -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto
  GetRightHandSideBuffer(InterfaceDescription const interface_description) const
      -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

  auto GetReinitializedBuffer(InterfaceDescription const interface_description)
      -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto
  GetReinitializedBuffer(InterfaceDescription const interface_description) const
      -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

  auto GetInitialBuffer(InterfaceDescription const interface_description)
      -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto GetInitialBuffer(InterfaceDescription const interface_description) const
      -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

  auto GetIntegratedBuffer(InterfaceDescription const interface_description)
      -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto
  GetIntegratedBuffer(InterfaceDescription const interface_description) const
      -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

  template <InterfaceDescriptionBufferType C>
  InterfaceDescriptions &GetInterfaceDescriptionBuffer();

  template <InterfaceDescriptionBufferType C>
  InterfaceDescriptions const &GetInterfaceDescriptionBuffer() const;

  InterfaceDescriptions &GetBaseBuffer();
  InterfaceDescriptions const &GetBaseBuffer() const;
  InterfaceDescriptions &GetRightHandSideBuffer();
  InterfaceDescriptions const &GetRightHandSideBuffer() const;
  InterfaceDescriptions &GetReinitializedBuffer();
  InterfaceDescriptions const &GetReinitializedBuffer() const;
  InterfaceDescriptions &GetInitialBuffer();
  InterfaceDescriptions const &GetInitialBuffer() const;
  InterfaceDescriptions &GetIntegratedBuffer();
  InterfaceDescriptions const &GetIntegratedBuffer() const;
  InterfaceDescriptions &GetInterfaceDescriptionBuffer(
      InterfaceDescriptionBufferType const buffer_type);
  InterfaceDescriptions const &GetInterfaceDescriptionBuffer(
      InterfaceDescriptionBufferType const buffer_type) const;

  // Returning state buffers
  auto GetInterfaceStateBuffer(InterfaceState const state_type)
      -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto GetInterfaceStateBuffer(InterfaceState const state_type) const
      -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

  InterfaceStates &GetInterfaceStateBuffer();
  InterfaceStates const &GetInterfaceStateBuffer() const;

  // returning parameter buffers
  auto GetInterfaceParameterBuffer(InterfaceParameter const parameter_type)
      -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto
  GetInterfaceParameterBuffer(InterfaceParameter const parameter_type) const
      -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

  InterfaceParameters &GetInterfaceParameterBuffer();
  InterfaceParameters const &GetInterfaceParameterBuffer() const;

  // returning general interface block buffer
  auto GetBuffer(InterfaceBlockBufferType const buffer_type)
      -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
  auto GetBuffer(InterfaceBlockBufferType const buffer_type) const
      -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
};

#endif // INTERFACE_BLOCK_H
