//===---------------------- interface_extender.h --------------------------===//
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
#ifndef INTERFACE_EXTENDER_H
#define INTERFACE_EXTENDER_H

#include "halo_manager.h"
#include "levelset/geometry/geometry_calculator_marching_cubes.h"
#include "user_specifications/numerical_setup.h"

/**
 * @brief The InterfaceExtender class extends interface fields in the narrow
 * band around the interface by solving the extension equation.
 *
 * @tparam DerivedInterfaceExtender Static Derived extender class which performs
 * the actual iterative extension
 * @tparam field_type The Interface field type for which the extender is used
 * (States or Parameters)
 */
template <typename DerivedInterfaceExtender, InterfaceFieldType field_type>
class InterfaceExtender {

  friend DerivedInterfaceExtender;

  // Class obtained from main
  HaloManager &halo_manager_;
  LogWriter &logger_;

  // private variables required for definition on which material field type the
  // extender works (static required for array definition)
  static constexpr InterfaceFieldType field_type_ = field_type;
  /**
   * The number of quantities that are necessary to track convergence. Those are
   * the maximum values of the quantities to extend and the residuum.
   */
  static constexpr unsigned int number_of_convergence_tracking_quantities_ =
      IF::NOFTE(field_type_) + 1;
  // numerical threshold
  static constexpr double epsilon_ = std::numeric_limits<double>::epsilon();

  /**
   * @brief Determines the maximum interface fields in the cells where we
   * extend. This is necessary to have a global normalization constant for
   * convergence tracking.
   * @param node The node for which the maximum interface fields are determined.
   * @param convergence_tracking_quantities A vector holding information about
   * the convergence status of the iterative extension method.
   */
  void DetermineMaximumValueOfQuantitiesToExtend(
      Node const &node,
      std::vector<double> &convergence_tracking_quantities) const {

    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

    for (unsigned int field_index = 0; field_index < IF::NOFTE(field_type_);
         ++field_index) {
      double const(&interface_field)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetFieldBuffer(
              field_type_, IF::FITE(field_type_, field_index));
      for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
        for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
          for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            if (std::abs(interface_tags[i][j][k]) < ITTI(IT::BulkPhase) &&
                std::abs(interface_tags[i][j][k]) > ITTI(IT::NewCutCell)) {
              convergence_tracking_quantities[field_index] =
                  std::max(convergence_tracking_quantities[field_index],
                           std::abs(interface_field[i][j][k]));
            }
          } // k
        }   // j
      }     // i
    }       // field of interface field_type
  }

  /**
   * @brief Default constructor of the InterfaceExtender class.
   * @param halo_manager Instance to a HaloManager which provides MPI-related
   * methods.
   */
  explicit InterfaceExtender(HaloManager &halo_manager)
      : halo_manager_(halo_manager), logger_(LogWriter::Instance()) {
    // Empty besides initializer list.
  }

public:
  InterfaceExtender() = delete;
  ~InterfaceExtender() = default;
  InterfaceExtender(InterfaceExtender const &) = delete;
  InterfaceExtender &operator=(InterfaceExtender const &) = delete;
  InterfaceExtender(InterfaceExtender &&) = delete;
  InterfaceExtender &operator=(InterfaceExtender &&) = delete;

  /**
   * @brief Performs an extension of interface quantities.
   * @param nodes The nodes for which scale separation should be done.
   */
  void Extend(std::vector<std::reference_wrapper<Node>> const &nodes) const {
    // Initialization of tracking quantities
    std::vector<double> convergence_tracking_quantities(
        number_of_convergence_tracking_quantities_, 0.0);
    // actual iterative loop
    for (unsigned int iteration_number = 0;
         iteration_number <
         InterfaceStateExtensionConstants::MaximumNumberOfIterations;
         ++iteration_number) {
      if constexpr (InterfaceStateExtensionConstants::TrackConvergence) {
        for (unsigned int field_index = 0; field_index < IF::NOFTE(field_type_);
             ++field_index) {
          convergence_tracking_quantities[field_index] = 0.0;
        }
        for (auto const &node : nodes) {
          DetermineMaximumValueOfQuantitiesToExtend(
              node, convergence_tracking_quantities);
        }

        MPI_Allreduce(MPI_IN_PLACE, convergence_tracking_quantities.data(),
                      number_of_convergence_tracking_quantities_, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);

        if (convergence_tracking_quantities[IF::NOFTE(field_type_)] <
                InterfaceStateExtensionConstants::MaximumResiduum &&
            iteration_number != 0) {
          if constexpr (GeneralTwoPhaseSettings::LogConvergenceInformation) {
            logger_.BufferMessage(
                "IntExt: " +
                std::to_string(static_cast<int>(iteration_number)) + " ");
          }
          break;
        } else if (iteration_number ==
                   InterfaceStateExtensionConstants::MaximumNumberOfIterations -
                       1) {
          if constexpr (GeneralTwoPhaseSettings::LogConvergenceInformation) {
            logger_.BufferMessage("IntExt: nc   !!!   ");
          }
        }
        convergence_tracking_quantities[IF::NOFTE(field_type_)] = 0.0;
      }

      // carry out the actual iterative extension on all nodes
      for (auto const &node : nodes) {
        static_cast<DerivedInterfaceExtender const &>(*this).IterativeExtension(
            node, convergence_tracking_quantities);
      }

      // Halo Update for all interface fields that should be extended
      for (unsigned int field_index = 0; field_index < IF::NOFTE(field_type_);
           ++field_index) {
        halo_manager_.InterfaceHaloUpdateOnLmax(
            MapInterfaceFieldToInterfaceBlockBufferType(
                field_type_, IF::FITE(field_type_, field_index)));
      }
    }
  }
};

#endif // INTERFACE_EXTENDER_H
