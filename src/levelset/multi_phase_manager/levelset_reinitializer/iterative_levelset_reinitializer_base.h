//===------------ iterative_levelset_reinitializer_base.h -----------------===//
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
#ifndef ITERATIVE_LEVELSET_REINITIALIZER_BASE_H
#define ITERATIVE_LEVELSET_REINITIALIZER_BASE_H

#include "enums/interface_tag_definition.h"
#include "levelset_reinitializer.h"
#include "user_specifications/two_phase_constants.h"
#include "utilities/buffer_operations_interface.h"
#include "utilities/mathematical_functions.h"

/**
 * @brief The class IterativeLevelsetReinitializerBase ensures the
 * (signed-)distance property of a level-set field.
 * @tparam Typename as template parameter due to CRTP.
 */
template <typename DerivedIterativeLevelsetReinitializer>
class IterativeLevelsetReinitializerBase
    : public LevelsetReinitializer<DerivedIterativeLevelsetReinitializer> {

  friend LevelsetReinitializer<DerivedIterativeLevelsetReinitializer>;
  friend DerivedIterativeLevelsetReinitializer;

  using LevelsetReinitializer<
      DerivedIterativeLevelsetReinitializer>::halo_manager_;
  using LevelsetReinitializer<DerivedIterativeLevelsetReinitializer>::logger_;

  /**
   * @brief The default constructor for a LevelsetReinitializer object.
   * @param halo_manager Instance to a HaloManager which provides MPI-related
   * methods.
   */
  explicit IterativeLevelsetReinitializerBase(HaloManager &halo_manager)
      : LevelsetReinitializer<DerivedIterativeLevelsetReinitializer>(
            halo_manager) {
    // Empty Constructor, besides initializer list.
  }

  /**
   * @brief Sets the cut-off in the levelset field of a node.
   * @param node The node with levelset block which has to be reinitialized.
   * @param levelset_type Level set buffer type which is reinitialized.
   */
  void
  CutOffSingleNode(Node &node,
                   InterfaceDescriptionBufferType const levelset_type) const {

    InterfaceBlock &interface_block = node.GetInterfaceBlock();
    double(&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        interface_block.GetInterfaceDescriptionBuffer(
            levelset_type)[InterfaceDescription::Levelset];

    // Cells which have a levelset value greater than the cutoff value are set
    // to cutoff.
    double const cutoff = CC::LSCOF();
    for (unsigned int i = 0; i < CC::TCX(); ++i) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int k = 0; k < CC::TCZ(); ++k) {
          if (std::abs(levelset[i][j][k]) > cutoff) {
            levelset[i][j][k] = Signum(levelset[i][j][k]) * cutoff;
          }
        } // k
      }   // j
    }     // i
  }

  /**
   * @brief Reinitializes a single-level set field as described in \cite
   * Sussman1994.
   * @param nodes Vector holding all nodes that have to be updated.
   * @param levelset_type Level set buffer type which is reinitialized.
   * @param is_last_stage whether it's the last RK stage or not.
   */
  void ReinitializeImplementation(
      std::vector<std::reference_wrapper<Node>> const &nodes,
      InterfaceDescriptionBufferType const levelset_type,
      bool const is_last_stage) const {

    InterfaceBlockBufferType const levelset_buffer_type =
        levelset_type == InterfaceDescriptionBufferType::Reinitialized
            ? InterfaceBlockBufferType::LevelsetReinitialized
            : InterfaceBlockBufferType::LevelsetIntegrated;
    // Store the original levelset field in the right-hand side buffer to have a
    // reference during reinitialization
    if (levelset_type == InterfaceDescriptionBufferType::Reinitialized) {
      BO::Interface::CopyInterfaceDescriptionBufferForNodeList<
          InterfaceDescriptionBufferType::Reinitialized,
          InterfaceDescriptionBufferType::RightHandSide>(nodes);
    } else {
      BO::Interface::CopyInterfaceDescriptionBufferForNodeList<
          InterfaceDescriptionBufferType::Integrated,
          InterfaceDescriptionBufferType::RightHandSide>(nodes);
    }

    // Carry out the actual reinitialization procedure
    double residuum = 0.0;
    for (unsigned int iteration_number = 0;
         iteration_number <
         ReinitializationConstants::MaximumNumberOfIterations;
         ++iteration_number) {

      if constexpr (ReinitializationConstants::TrackConvergence) {
        residuum = 0.0;
      }
      for (auto &node : nodes) {
        residuum = std::max(
            residuum,
            static_cast<DerivedIterativeLevelsetReinitializer const &>(*this)
                .ReinitializeSingleNodeImplementation(node, levelset_type,
                                                      is_last_stage));
      }

      // halo update
      halo_manager_.InterfaceHaloUpdateOnLmax(levelset_buffer_type);
      if constexpr (ReinitializationConstants::TrackConvergence) {
        MPI_Allreduce(MPI_IN_PLACE, &residuum, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);

        if (residuum < ReinitializationConstants::MaximumResiduum) {
          if constexpr (GeneralTwoPhaseSettings::LogConvergenceInformation) {
            logger_.BufferMessage(
                "Reinit: " +
                std::to_string(static_cast<int>(iteration_number)) + " ");
          }
          break;
        } else if (iteration_number ==
                   ReinitializationConstants::MaximumNumberOfIterations - 1) {
          if constexpr (GeneralTwoPhaseSettings::LogConvergenceInformation) {
            logger_.BufferMessage("Reinit: nc   !!!   ");
          }
        }
      }
    }

    for (auto &node : nodes) {
      CutOffSingleNode(node, levelset_type);
    }
    halo_manager_.InterfaceHaloUpdateOnLmax(levelset_buffer_type);
  }

public:
  IterativeLevelsetReinitializerBase() = delete;
  virtual ~IterativeLevelsetReinitializerBase() = default;
  IterativeLevelsetReinitializerBase(
      IterativeLevelsetReinitializerBase const &) = delete;
  IterativeLevelsetReinitializerBase &
  operator=(IterativeLevelsetReinitializerBase const &) = delete;
  IterativeLevelsetReinitializerBase(IterativeLevelsetReinitializerBase &&) =
      delete;
  IterativeLevelsetReinitializerBase &
  operator=(IterativeLevelsetReinitializerBase &&) = delete;
};

#endif // ITERATIVE_LEVELSET_REINITIALIZER_BASE_H
