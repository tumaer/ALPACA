//===------------------ two_phase_interface_extender.h --------------------===//
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
#ifndef TWO_PHASE_INTERFACE_EXTENDER_H
#define TWO_PHASE_INTERFACE_EXTENDER_H

#include "interface_extender.h"

template <InterfaceFieldType field_type>
class TwoPhaseInterfaceExtender
    : public InterfaceExtender<TwoPhaseInterfaceExtender<field_type>,
                               field_type> {
  // Definition for shortening
  using InterfaceExtenderSpecification =
      InterfaceExtender<TwoPhaseInterfaceExtender<field_type>, field_type>;

  friend InterfaceExtenderSpecification;

private:
  // use some variables from base class
  using InterfaceExtenderSpecification::epsilon_;
  using InterfaceExtenderSpecification::field_type_;
  using InterfaceExtenderSpecification::
      number_of_convergence_tracking_quantities_;
  static constexpr unsigned int stencil_width_ = 2;
  static constexpr unsigned int i_offset_ = CC::HS() - stencil_width_;
  static constexpr unsigned int j_offset_ =
      CC::DIM() != Dimension::One ? CC::HS() - stencil_width_ : 0;
  static constexpr unsigned int k_offset_ =
      CC::DIM() == Dimension::Three ? CC::HS() - stencil_width_ : 0;
  static constexpr unsigned int repetition_ = CC::HS() / stencil_width_;

  /**
   * @brief Actual implementation of the iterative extension equation to extend
   * interface fields into narrow band.
   * @param node The node which contains the phase which is extended.
   * @param convergence_tracking_quantities A vector holding information about
   * the convergence status of the iterative extension method.
   */
  void IterativeExtension(
      Node &node, std::vector<double> &convergence_tracking_quantities) const {

    double const(&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceBlock().GetReinitializedBuffer(
            InterfaceDescription::Levelset);
    std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] =
        node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

    std::array<double, DTI(CC::DIM())> rhs_contributions;

    double interface_field_change[CC::TCX()][CC::TCY()][CC::TCZ()];

    for (unsigned int field_index = 0; field_index < IF::NOFTE(field_type_);
         ++field_index) {
      double const one_normalization_constant =
          1.0 / (convergence_tracking_quantities[field_index] + epsilon_);
      double(&interface_field)[CC::TCX()][CC::TCY()][CC::TCZ()] =
          node.GetInterfaceBlock().GetFieldBuffer(
              field_type_, IF::FITE(field_type_, field_index));

      for (unsigned int iteration = 0; iteration < repetition_; ++iteration) {
        for (unsigned int i = 0; i < CC::TCX(); ++i) {
          for (unsigned int j = 0; j < CC::TCY(); ++j) {
            for (unsigned int k = 0; k < CC::TCZ(); ++k) {
              interface_field_change[i][j][k] = 0.0;
            } // k
          }   // j
        }     // i

        // Loop through all narrow-band cells which are no cut-cells to compute
        // increment for iterative reinitialization
        for (unsigned int i = CC::FICX() - i_offset_;
             i <= CC::LICX() + i_offset_; ++i) {
          for (unsigned int j = CC::FICY() - j_offset_;
               j <= CC::LICY() + j_offset_; ++j) {
            for (unsigned int k = CC::FICZ() - k_offset_;
                 k <= CC::LICZ() + k_offset_; ++k) {
              /**
               * We also extend in the reinitialization band in order to have
               * better convergence behaviour (Reduce influence of implicitly
               * imposed boundary conditions at the end of the narrow band). The
               * convergence criteria is only checked for the extension band.
               */
              if (std::abs(interface_tags[i][j][k]) <=
                      ITTI(IT::ReinitializationBand) &&
                  std::abs(interface_tags[i][j][k]) > ITTI(IT::NewCutCell)) {

                // compute normal
                std::array<double, 3> const normal =
                    GetNormal(levelset_reinitialized, i, j, k);

                double const ls_sign = Signum(levelset_reinitialized[i][j][k]);

                // compute extension terms and store in vector
                // i-direction
                // direction of the gradient from interface towards the material
                // in both directions
                int const n_i =
                    (normal[0] * ls_sign < 0.0) - (normal[0] * ls_sign > 0.0);
                rhs_contributions[0] =
                    (interface_field[i + n_i][j][k] - interface_field[i][j][k]);
                rhs_contributions[0] *= std::abs(normal[0]);

                // j-direction
                if (CC::DIM() != Dimension::One) {
                  int const n_j =
                      (normal[1] * ls_sign < 0.0) - (normal[1] * ls_sign > 0.0);
                  rhs_contributions[1] = (interface_field[i][j + n_j][k] -
                                          interface_field[i][j][k]);
                  rhs_contributions[1] *= std::abs(normal[1]);
                }

                // k-direction
                if (CC::DIM() == Dimension::Three) {
                  int const n_k =
                      (normal[2] * ls_sign < 0.0) - (normal[2] * ls_sign > 0.0);
                  rhs_contributions[2] = (interface_field[i][j][k + n_k] -
                                          interface_field[i][j][k]);
                  rhs_contributions[2] *= std::abs(normal[2]);
                }

                interface_field_change[i][j][k] =
                    ConsistencyManagedSum(rhs_contributions) *
                    InterfaceStateExtensionConstants::Dtau;
                if (InterfaceStateExtensionConstants::TrackConvergence &&
                    std::abs(interface_tags[i][j][k]) <=
                        ITTI(IT::ExtensionBand)) {
                  convergence_tracking_quantities[IF::NOFTE(field_type_)] =
                      std::max(convergence_tracking_quantities[IF::NOFTE(
                                   field_type_)],
                               std::abs(interface_field_change[i][j][k] *
                                        one_normalization_constant));
                }
              } // if cut cell
            }   // k
          }     // j
        }       // i

        for (unsigned int i = CC::FICX() - i_offset_;
             i <= CC::LICX() + i_offset_; ++i) {
          for (unsigned int j = CC::FICY() - j_offset_;
               j <= CC::LICY() + j_offset_; ++j) {
            for (unsigned int k = CC::FICZ() - k_offset_;
                 k <= CC::LICZ() + k_offset_; ++k) {
              interface_field[i][j][k] += interface_field_change[i][j][k];
            } // k
          }   // j
        }     // i
      }

    } // field of interface field_type
  }

public:
  TwoPhaseInterfaceExtender() = delete;
  explicit TwoPhaseInterfaceExtender(HaloManager &halo_manager)
      : InterfaceExtenderSpecification(halo_manager) {
    /** Empty besides initilializer list for call of base class constructor */
  }
  ~TwoPhaseInterfaceExtender() = default;
  TwoPhaseInterfaceExtender(TwoPhaseInterfaceExtender const &) = delete;
  TwoPhaseInterfaceExtender &
  operator=(TwoPhaseInterfaceExtender const &) = delete;
  TwoPhaseInterfaceExtender(TwoPhaseInterfaceExtender &&) = delete;
  TwoPhaseInterfaceExtender &operator=(TwoPhaseInterfaceExtender &&) = delete;
};

#endif // TWO_PHASE_INTERFACE_EXTENDER_H
