//===---------------------- ghost_fluid_extender.h ------------------------===//
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
#ifndef GHOST_FLUID_EXTENDER_H
#define GHOST_FLUID_EXTENDER_H

#include "halo_manager.h"
#include "levelset/geometry/geometry_calculator_marching_cubes.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "materials/material_manager.h"
#include "user_specifications/numerical_setup.h"
#include "user_specifications/two_phase_constants.h"

/**
 * @brief The GhostFluidExtender class extends material-states from
 * real-material to ghost-material by iteratively solving the extension
 * equation.
 *
 * @tparam DerivedGhostFluidExtender Static Derived extender class which
 * performs the actual iterative extension
 * @tparam field_type The Fluid field type for which the extender is used
 * (PrimeStates, Parameters or Conservatives)
 */
template <typename DerivedGhostFluidExtender, MaterialFieldType field_type>
class GhostFluidExtender {

  friend DerivedGhostFluidExtender;

private:
  // member classes
  MaterialManager const &material_manager_;
  HaloManager
      &halo_manager_; // TODO-19 NH Think about making it const (rats tail)
  LogWriter &logger_;

  // private variables required for definition on which material field type the
  // extender works (static required for array definition)
  static constexpr MaterialFieldType field_type_ = field_type;
  /**
   * The number of quantities that are necessary to track convergence. Those are
   * the maximum values of the quantities to extend and the residuum.
   */
  static constexpr unsigned int number_of_convergence_tracking_quantities_ =
      MF::ANOF(field_type_) + 1;
  // numerical threshold
  static constexpr double epsilon_ = std::numeric_limits<double>::epsilon();

  /**
   * @brief Determines the maximum material fields in the cells where we extend.
   * This is necessary to have a global normalization constant for convergence
   * tracking.
   * @param node The node for which the maximum material fields are determined.
   * @param convergence_tracking_quantities An array holding information about
   * the convergence status of the iterative extension method.
   */
  void DetermineMaximumValueOfQuantitiesToExtend(
      Node const &node,
      double (&convergence_tracking_quantities)
          [2][number_of_convergence_tracking_quantities_]) const {

    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();
    double const(&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetReinitializedBuffer(
            InterfaceDescription::VolumeFraction);

    // Loop through all materials on the node
    for (auto const &phase : node.GetPhases()) {
      // Define the material specific parameter
      std::int8_t const material_sign =
          MaterialSignCapsule::SignOfMaterial(phase.first);
      unsigned int const material_index =
          phase.first == MaterialSignCapsule::PositiveMaterial() ? 0 : 1;
      double const reference_volume_fraction = (material_sign > 0) ? 0.0 : 1.0;
      double const material_sign_double = double(material_sign);
      // Loop through all fields of the given field_Type
      for (unsigned int field_index = 0; field_index < MF::ANOF(field_type_);
           field_index++) {
        double const(&extension_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] =
            phase.second.GetFieldBuffer(field_type_, field_index);
        // Loop through all internal cells of the block
        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
          for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
            for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
              double const cell_volume_fraction =
                  reference_volume_fraction +
                  material_sign_double * volume_fraction[i][j][k];
              // Change the convergence of each field
              if (std::abs(interface_tags[i][j][k]) <=
                      ITTI(IT::ReinitializationBand) &&
                  (cell_volume_fraction <= CC::ETH())) {
                convergence_tracking_quantities[material_index][field_index] =
                    std::max(convergence_tracking_quantities[material_index]
                                                            [field_index],
                             std::abs(extension_buffer[i][j][k]));
              }
            } // k
          }   // j
        }     // i
      }       // extension_buffer
    }         // phases
  }

  /**
   * @brief The default constructor for the IterativeGhostFluidExtender. Calls
   * the default constructor of the base class.
   * @param material_manager Instance of a material manager, which already has
   * been initialized according to the user input.
   * @param halo_manager Instance to a HaloManager which provides MPI-related
   * methods.
   */
  explicit GhostFluidExtender(MaterialManager const &material_manager,
                              HaloManager &halo_manager)
      : material_manager_(material_manager), halo_manager_(halo_manager),
        logger_(LogWriter::Instance()) {
    /** Empty besides initializer list */
  }

public:
  GhostFluidExtender() = delete;
  ~GhostFluidExtender() = default;
  GhostFluidExtender(GhostFluidExtender const &) = delete;
  GhostFluidExtender &operator=(GhostFluidExtender const &) = delete;
  GhostFluidExtender(GhostFluidExtender &&) = delete;
  GhostFluidExtender operator=(GhostFluidExtender &&) = delete;

  /**
   * @brief Iteratively solves the extension equation. For description of
   * functionality also see base class.
   * @param nodes The nodes for which the extension equation is solved.
   */
  void Extend(std::vector<std::reference_wrapper<Node>> const &nodes) const {

    // The 2 is hardcoded on purpose. It corresponds to the number of materials.
    // An issue about that is already in the git.
    double convergence_tracking_quantities
        [2][number_of_convergence_tracking_quantities_];
    // Initialize the array with zeros
    for (unsigned int material_index = 0; material_index < 2;
         ++material_index) {
      for (unsigned int field_index = 0; field_index <= MF::ANOF(field_type_);
           ++field_index) {
        convergence_tracking_quantities[material_index][field_index] = 0.0;
      }
    }

    // Actual iterative loop
    for (unsigned int iteration_number = 0;
         iteration_number < ExtensionConstants::MaximumNumberOfIterations;
         ++iteration_number) {
      // Additional computation if the convergence is tracked
      if constexpr (ExtensionConstants::TrackConvergence) {
        // Reset tracking quantities
        for (unsigned int material_index = 0; material_index < 2;
             ++material_index) {
          for (unsigned int field_index = 0;
               field_index < MF::ANOF(field_type_); ++field_index) {
            convergence_tracking_quantities[material_index][field_index] = 0.0;
          }
        }
        for (Node const &node : nodes) {
          DetermineMaximumValueOfQuantitiesToExtend(
              node, convergence_tracking_quantities);
        }

        MPI_Allreduce(MPI_IN_PLACE, &convergence_tracking_quantities,
                      number_of_convergence_tracking_quantities_ * 2,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        // Write convergence to logger if desired or if maximum of iterations is
        // reached
        if (convergence_tracking_quantities[0][MF::ANOF(field_type_)] <
                ExtensionConstants::MaximumResiduum &&
            convergence_tracking_quantities[1][MF::ANOF(field_type_)] <
                ExtensionConstants::MaximumResiduum &&
            iteration_number != 0) {
          if constexpr (GeneralTwoPhaseSettings::LogConvergenceInformation) {
            logger_.BufferMessage(
                "Ext: " + std::to_string(static_cast<int>(iteration_number)) +
                " ");
          }
          break;
        } else if (iteration_number ==
                   ExtensionConstants::MaximumNumberOfIterations - 1) {
          if constexpr (GeneralTwoPhaseSettings::LogConvergenceInformation) {
            logger_.BufferMessage("Ext: nc   !!!   ");
          }
        }
        // Reset the maximum tracking quantity
        convergence_tracking_quantities[0][MF::ANOF(field_type_)] = 0.0;
        convergence_tracking_quantities[1][MF::ANOF(field_type_)] = 0.0;
      }

      // iterative extension on field buffer (static derived extender)
      for (Node &node : nodes) {
        static_cast<DerivedGhostFluidExtender const &>(*this)
            .IterativeExtension(node, convergence_tracking_quantities);
      } // nodes

      // Update the halos after each iterative step
      halo_manager_.MaterialHaloUpdateOnLmaxMultis(field_type_);
    }
  };
};

#endif // GHOST_FLUID_EXTENDER_H
