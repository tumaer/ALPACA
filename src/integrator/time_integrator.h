//===----------------------- time_integrator.h ----------------------------===//
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
#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H

#include <numeric>

#include "block_definitions/block.h"
#include "boundary_condition/boundary_specifications.h"
#include "enums/interface_tag_definition.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "topology/node.h"
#include "utilities/buffer_operations_interface.h"
#include "utilities/buffer_operations_material.h"

/**
 * @brief The TimeIntegrator class defines an interface for the time integration
 * of the underlying equations. TimeIntegrator (sub) classes provide
 * scheme-specific parameters.
 * @tparam DerivedTimeIntegrator The template for the derived classes. This is
 * necessary for the CRTP.
 */
template <typename DerivedTimeIntegrator> class TimeIntegrator {

  friend DerivedTimeIntegrator;

  double current_run_time_;
  std::vector<double> micro_timestep_sizes_;

  /**
   * @brief Constructor.
   * @param start_time Time when the simulation should start.
   */
  explicit TimeIntegrator(double const start_time = 0.0)
      : current_run_time_(start_time) {}

  /**
   * @brief Performs the time integration (same Algorithm as in Integrate
   * function) for the Jump Flux Buffers.
   * @param block The block whose contents should be integrated, consistent
   * buffer states need to be ensured by caller.
   * @param timestep The size of the time step used in the current integration
   * step.
   */
  void IntegrateJumpConservatives(Block &block, double const timestep) const {

    for (auto const &location : CC::ANBS()) {
      double(&boundary_conservatives)[MF::ANOE()][CC::ICY()][CC::ICZ()] =
          block.GetBoundaryJumpConservatives(location);
      double(&boundary_fluxes)[MF::ANOE()][CC::ICY()][CC::ICZ()] =
          block.GetBoundaryJumpFluxes(location);
      for (unsigned int e = 0; e < MF::ANOE(); ++e) {
        for (unsigned int i = 0; i < CC::ICY(); ++i) {
          for (unsigned int j = 0; j < CC::ICZ(); ++j) {
            // integrate change of conservatives over the block boundary
            boundary_conservatives[e][i][j] +=
                timestep * boundary_fluxes[e][i][j];

            // reset boundary fluxes to 0
            boundary_fluxes[e][i][j] = 0.0;
          }
        }
      }
    }
  }

  /**
   * @brief Increments the current solution by one timestep (or stage for
   * Runge-Kutta methods). Does not perform correctness checks, the provided
   * block buffers must be in a consistent state with the respective integration
   * scheme.
   * @param block The block whose contents should be integrated, consistent
   * buffer states need to be ensured by caller.
   * @param timestep The size of the time step used in the current integration
   * step.
   * @note  This function works only for RK2 and RK3. In order to use
   * higher-order schemes, it might be necessary to adapt it.
   */
  void IntegrateConservatives(Block &block, double const timestep) const {

    for (Equation const eq : MF::ASOE()) {
      double(&u_old)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          block.GetAverageBuffer(eq);
      double(&u_new)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          block.GetRightHandSideBuffer(eq);
      for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
        for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
          for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            u_new[i][j][k] = u_old[i][j][k] + timestep * u_new[i][j][k];
          }
        }
      }
    }
  }

  /**
   * @brief Same as for IntegrateConservatives(Block& block, double const
   * timestep) function, however, solution is only incremented in time in the
   * halo cells. Integrates all (i.e. six in three dimensions) halo cells.  Does
   * not perform correctness checks, the provided block buffers must be in a
   * consistent state with the respective integration scheme.
   * @param block The block whose contents should be integrated, consistent
   * buffer states need to be ensured by caller.
   * @param timestep The size of the time step used in the current integration
   * step.
   * @param start_indices_halo The start indices of the halo cells in all
   * directions.
   * @param halo_size The number of halo cells to be integrated in each
   * direction.
   */
  void IntegrateHalo(Block &block, double const timestep,
                     std::array<int, 3> const start_indices_halo,
                     std::array<int, 3> const halo_size) const {

    for (Equation const eq : MF::ASOE()) {
      double(&u_old)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          block.GetAverageBuffer(eq);
      double(&u_new)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          block.GetRightHandSideBuffer(eq);

      for (int i = start_indices_halo[0];
           i < start_indices_halo[0] + halo_size[0]; i++) {
        for (int j = start_indices_halo[1];
             j < start_indices_halo[1] + halo_size[1]; j++) {
          for (int k = start_indices_halo[2];
               k < start_indices_halo[2] + halo_size[2]; k++) {
            u_new[i][j][k] = u_old[i][j][k] + timestep * u_new[i][j][k];
          }
        }
      }
    }
  }

  /**
   * @brief Increments the current level-set field by one timestep (or stage for
   * Runge-Kutta methods). Does not perform correctness checks, the provided
   * block buffers must be in a consistent state with the respective integration
   * scheme.
   * @param node The node whose level-set field should be incremented.
   * @param timestep The size of the time step used in the current integration
   * step.
   * @note  This function works only for RK2 and RK3. In order to use
   * higher-order schemes, it might be necessary to adapt it.
   */
  void IntegrateLevelset(Node &node, double const timestep) const {
    double(&levelset_new)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetRightHandSideBuffer(
            InterfaceDescription::Levelset);
    double const(&levelset_old)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetBaseBuffer(InterfaceDescription::Levelset);
    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
    double const(&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetReinitializedBuffer(
            InterfaceDescription::Levelset);

    for (unsigned int i = 0; i < CC::TCX(); ++i) {
      for (unsigned int j = 0; j < CC::TCY(); ++j) {
        for (unsigned int k = 0; k < CC::TCZ(); ++k) {
          // integrate narrow band including cut-cells
          if (std::abs(interface_tags[i][j][k]) < ITTI(IT::BulkPhase)) {
            levelset_new[i][j][k] =
                levelset_old[i][j][k] + timestep * levelset_new[i][j][k];
          } else {
            levelset_new[i][j][k] = levelset_reinitialized[i][j][k];
          }
        } // k
      }   // j
    }     // i
  }

  /**
   * @brief Returns the timestep size multiplication factor of the jump
   * conservatives for the current RK stage.
   * @param stage The current stage of the RK scheme.
   */
  constexpr double
  GetTimestepMultiplierJumpConservatives(unsigned int const stage) const {
    return DerivedTimeIntegrator::timestep_multiplier_jump_conservatives_
        [stage];
  }

  /**
   * @brief Returns the timestep size multiplication factor of the conservatives
   * for the current RK stage.
   * @param stage The current stage of the RK scheme.
   */
  constexpr double
  GetTimestepMultiplierConservatives(unsigned int const stage) const {
    return DerivedTimeIntegrator::timestep_multiplier_conservatives_[stage];
  }

  /**
   * @brief Returns the buffer multiplication factor of the conservatives for
   * the current RK stage.
   * @param stage The current stage of the RK scheme.
   */
  constexpr auto GetBufferMultiplier(unsigned int const stage) const {
    return DerivedTimeIntegrator::buffer_multiplier_[stage - 1];
  }

public:
  TimeIntegrator() = delete;
  ~TimeIntegrator() = default;
  TimeIntegrator(TimeIntegrator const &) = delete;
  TimeIntegrator &operator=(TimeIntegrator const &) = delete;
  TimeIntegrator(TimeIntegrator &&) = delete;
  TimeIntegrator &operator=(TimeIntegrator &&) = delete;

  /**
   * @brief Sets the start time of the simulation to the given value.
   * @param The start time to use for the simlation.
   */
  void SetStartTime(double const start_time) { current_run_time_ = start_time; }

  /**
   * @brief Adds a micro time step (size) to the list of timestep_sizes on
   * finest level. $NOT SAFE: Corrupted Input results in wrong integrations$
   * @param time_step The time step (size) to be appended.
   */
  void AppendMicroTimestep(double const time_step) {
    micro_timestep_sizes_.push_back(time_step);
  }

  /**
   * @brief Computes the macro time step size, adds it to the macro timestep
   * list and empties the micro timestep list.
   */
  void FinishMacroTimestep() {
    current_run_time_ += std::accumulate(micro_timestep_sizes_.cbegin(),
                                         micro_timestep_sizes_.cend(), 0.0);
    micro_timestep_sizes_.clear();
  }

  /**
   * @brief Gives the current list of micro time step sizes.
   * @return time step sizes on the finest level.
   */
  std::vector<double> const &MicroTimestepSizes() const {
    return micro_timestep_sizes_;
  }

  /**
   * @brief Returns the current run time, i.e. time of all fully passed MACRO
   * timesteps.
   * @return Run time.
   */
  inline double CurrentRunTime() const { return current_run_time_; }

  /**
   * @brief Integrates all jump halos a node holds.
   * @param node The node under consideration, may not have jumps.
   * @param stage The integration stage.
   * @param number_of_timesteps The number of time steps relevant for this
   * integration, i.e. on coarser levels the timestep sizes of the finer levels
   * need to be summed.
   * @param start_indices_halo The start indices of the halo cells in all
   * directions.
   * @param halo_size The number of halo cells to be integrated in each
   * direction.
   */
  void IntegrateJumpHalos(Node &node, unsigned int const stage,
                          unsigned int const number_of_timesteps,
                          std::array<int, 3> const start_indices_halo,
                          std::array<int, 3> const halo_size) const {
#ifndef PERFORMANCE
    if (stage >= NumberOfStages()) {
      throw std::invalid_argument(
          "Stage is too large for the chosen time integration scheme");
    }
#endif

    double const timestep =
        std::accumulate(micro_timestep_sizes_.crbegin(),
                        micro_timestep_sizes_.crbegin() + number_of_timesteps,
                        0.0) *
        GetTimestepMultiplierConservatives(stage);

    for (auto &phase : node.GetPhases()) {
      IntegrateHalo(phase.second, timestep, start_indices_halo, halo_size);
    }
  }

  /**
   * @brief Integrates one node by one stage. Integration is done for the
   * conservatives.
   * @param node The node to be integrated.
   * @param stage The integration stage.
   * @param number_of_timesteps The number of time steps relevant for this
   * integration, i.e. on coarser levels the timestep sizes of the finer levels
   * need to be summed.
   */
  void IntegrateNode(Node &node, unsigned int const stage,
                     unsigned int const number_of_timesteps) const {
#ifndef PERFORMANCE
    if (stage >= NumberOfStages()) {
      throw std::invalid_argument(
          "Stage is too large for the chosen time integration scheme");
    }
#endif

    // The timestepsize for the jump conservatives is different than for the
    // regular RK stage, depending on the order of the overall scheme.
    // Therefore, different multiplicators are required.
    double const timestep = std::accumulate(
        micro_timestep_sizes_.crbegin(),
        micro_timestep_sizes_.crbegin() + number_of_timesteps, 0.0);
    double const multiplier_jump_conservatives =
        GetTimestepMultiplierJumpConservatives(stage);
    double const multiplier_conservatives =
        GetTimestepMultiplierConservatives(stage);

    for (auto &phase : node.GetPhases()) {
      IntegrateJumpConservatives(phase.second,
                                 multiplier_jump_conservatives * timestep);
    }

    for (auto &phase : node.GetPhases()) {
      IntegrateConservatives(phase.second, multiplier_conservatives * timestep);
    }
  }

  /**
   * @brief Integrates one node by one stage. Integration is done for the
   * level-set field.
   * @param node The node to be integrated.
   * @param stage The integration stage.
   */
  void IntegrateLevelsetNode(Node &node, unsigned int const stage) const {
#ifndef PERFORMANCE
    if (stage >= NumberOfStages()) {
      throw std::invalid_argument(
          "Stage is too large for Runge-Kutta-3 integration");
    }
#endif
    double const timestep = micro_timestep_sizes_.back();
    double const multiplier_conservatives =
        GetTimestepMultiplierConservatives(stage);

    IntegrateLevelset(node, multiplier_conservatives * timestep);
  }

  /**
   * @brief Fills the initial buffer with values. This is necessary to realize
   * some time-stepping schemes.
   * @param node The node from which the necessary conservatives are written to
   * the initial buffer.
   * @param stage The stage of the time stepping scheme.
   */
  void FillInitialBuffer(Node &node, unsigned int const stage) const {
    if (stage == 0) {
      if (node.HasLevelset()) {
        for (auto &mat_block : node.GetPhases()) {
          std::int8_t const material_sign =
              MaterialSignCapsule::SignOfMaterial(mat_block.first);
          double const(&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] =
              node.GetInterfaceBlock().GetReinitializedBuffer(
                  InterfaceDescription::VolumeFraction);

          double const reference_volume_fraction =
              (material_sign > 0) ? 0.0 : 1.0;
          double const material_sign_double = double(material_sign);

          for (Equation const eq : MF::ASOE()) {
            double const(&u)[CC::TCX()][CC::TCY()][CC::TCZ()] =
                mat_block.second.GetAverageBuffer(eq);
            double(&u_initial)[CC::TCX()][CC::TCY()][CC::TCZ()] =
                mat_block.second.GetInitialBuffer(eq);
            for (unsigned int i = 0; i < CC::TCX(); ++i) {
              for (unsigned int j = 0; j < CC::TCY(); ++j) {
                for (unsigned int k = 0; k < CC::TCZ(); ++k) {
                  u_initial[i][j][k] =
                      u[i][j][k] *
                      (reference_volume_fraction +
                       material_sign_double * volume_fraction[i][j][k]);
                } // k
              }   // j
            }     // i
          }       // equations
        }         // phases
      } else {    // nodes without levelset
        BO::Material::CopyConservativeBuffersForNode<
            ConservativeBufferType::Average, ConservativeBufferType::Initial>(
            node);
      }
    } // if initial stage
  }

  /**
   * @brief Fills the initial level-set buffer with the initial values of the
   * RK-stage. This is only done for the first Runge-Kutta step.
   * @param node The node for which the level set of the last RK stage are
   * written to the temporary buffer.
   * @param stage The stage of the time stepping scheme.
   */
  void FillInitialLevelsetBuffer(Node &node, unsigned int const stage) const {
    if (stage == 0) {
      BO::Interface::CopyInterfaceDescriptionBufferForNode<
          InterfaceDescriptionBufferType::Base,
          InterfaceDescriptionBufferType::Initial,
          InterfaceDescription::Levelset>(node);
      BO::Interface::CopyInterfaceDescriptionBufferForNode<
          InterfaceDescriptionBufferType::Reinitialized,
          InterfaceDescriptionBufferType::Initial,
          InterfaceDescription::VolumeFraction>(node);
    } // stage
  }

  /**
   * @brief Prepares the average buffer for the next Runge-Kutta stage.
   * @param node The node for which the average buffer is prepared for the next
   * stage.
   * @param stage The stage of the time-stepping scheme.
   */
  void PrepareBufferForIntegration(Node &node, unsigned int const stage) const {
    // buffer preparation is only necessary for later stages
    if (stage != 0) {

      auto const multipliers = GetBufferMultiplier(stage);

      for (auto &mat_block : node.GetPhases()) {
        for (Equation const eq : MF::ASOE()) {
          double(&u)[CC::TCX()][CC::TCY()][CC::TCZ()] =
              mat_block.second.GetAverageBuffer(eq);
          double const(&u_initial)[CC::TCX()][CC::TCY()][CC::TCZ()] =
              mat_block.second.GetInitialBuffer(eq);
          for (unsigned int i = 0; i < CC::TCX(); ++i) {
            for (unsigned int j = 0; j < CC::TCY(); ++j) {
              for (unsigned int k = 0; k < CC::TCZ(); ++k) {
                u[i][j][k] = multipliers[0] * u[i][j][k] +
                             multipliers[1] * u_initial[i][j][k];
              } // k
            }   // j
          }     // i
        }       // equations
      }         // phases
    }
  }

  /**
   * @brief Prepares the level-set buffer for the next Runge-Kutta stage.
   * @param node The node for which the average buffer is prepared for the next
   * stage.
   * @param stage The stage of the time-stepping scheme.
   */
  void PrepareLevelsetBufferForIntegration(Node &node,
                                           unsigned int const stage) const {
    // buffer preparation is only necessary for intermediate and last stage
    if (stage != 0) {

      auto const multipliers = GetBufferMultiplier(stage);

      double(&levelset)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetBaseBuffer(
              InterfaceDescription::Levelset);
      double const(&levelset_initial)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetInitialBuffer(
              InterfaceDescription::Levelset);

      for (unsigned int i = 0; i < CC::TCX(); ++i) {
        for (unsigned int j = 0; j < CC::TCY(); ++j) {
          for (unsigned int k = 0; k < CC::TCZ(); ++k) {
            levelset[i][j][k] = multipliers[0] * levelset[i][j][k] +
                                multipliers[1] * levelset_initial[i][j][k];
          } // k
        }   // j
      }     // i
    }
  }

  /**
   * @brief Gives the number of stages performed in one timestep for this kind
   * of solver
   * @return Number of stages.
   */
  constexpr unsigned int NumberOfStages() const {
    return DerivedTimeIntegrator::number_of_stages_;
  }

  /**
   * @brief Gives whether or not the given stage is the last stage of this time
   * integrator.
   * @param stage The stage to be checked (zero based).
   * @return True if given stage is the last stage, false otherwise.
   */
  constexpr bool IsLastStage(unsigned int const stage) const {
    return stage == (NumberOfStages() - 1);
  }

  /**
   * @brief Swaps the two buffers of the provided Block object.
   * @param block Reference to Block object whose buffers are to be swapped.
   * @note This function must be inherited properly as soon as other integrators
   * are implemented!
   */
  void SwapBuffersForNextStage(Node &node) const {
    // swap the conservative buffers
    BO::Material::SwapConservativeBuffersForNode<
        ConservativeBufferType::RightHandSide, ConservativeBufferType::Average>(
        node);
    // swap the levelset buffer properly if required
    if (node.HasLevelset()) {
      BO::Interface::SwapInterfaceDescriptionBufferForNode<
          InterfaceDescriptionBufferType::RightHandSide,
          InterfaceDescriptionBufferType::Base>(node);
    }
  }
};

#endif // TIME_INTEGRATOR_H
