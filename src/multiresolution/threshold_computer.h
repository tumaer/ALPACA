//===---------------------- threshold_computer.h --------------------------===//
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
#ifndef THRESHOLD_COMPUTER_H
#define THRESHOLD_COMPUTER_H

#include "enums/dimension_definition.h"
#include "user_specifications/compile_time_constants.h"
#include <cmath>

/**
 * @brief A class to provide the multiresolution threshold values for the given
 * intput. According to \cite Harten1995 and \cite Roussel2003.
 * @tparam DIM Dimension of the current simulation.
 * @tparam A Limiting convergence order of the overall scheme (space + time).
 */
template <Dimension DIM, unsigned int A> class ThresholdComputer {

  static constexpr double D_ =
      static_cast<double>(DTI(DIM)); // Capitalization according to paper
  static constexpr double alpha_ = static_cast<double>(A);
  unsigned int const maximum_level_;
  double const reference_epsilon_;

  /**
   * @brief Computes the reference epsilon used in the multiresolution analysis
   * according to \cite Roussel2013.
   * @param user_epsilon_reference The threshold value the user wants to ensure
   * on the given level.
   * @param user_level_of_reference The level at which the given threshold
   * should be ensured.
   * @return Epsilon Reference to be used in multiresolution analysis.
   */
  inline double
  ComputeReferenceEpsilon(double const user_epsilon_reference,
                          unsigned int const user_level_of_reference) const {
    return user_epsilon_reference *
           std::pow(2.0, -1.0 * (alpha_ + 1.0) *
                             (double(maximum_level_) -
                              double(user_level_of_reference)));
  }

public:
  /**
   * @brief Standard constructor.
   * @param maximum_level Maximum level in the simulation.
   * @param user_reference_level Reference level provided by user input for the
   * multiresolution threshold.
   * @param user_reference_epsilon Reference epsilon provided by user input for
   * the multiresolution threshold.
   */
  explicit ThresholdComputer(unsigned int const maximum_level,
                             unsigned int const user_reference_level,
                             double const user_reference_epsilon)
      : maximum_level_(maximum_level),
        reference_epsilon_(ComputeReferenceEpsilon(user_reference_epsilon,
                                                   user_reference_level)) {
    // Empty besides initializer list
  }
  ThresholdComputer() = delete;
  ~ThresholdComputer() = default;
  ThresholdComputer(ThresholdComputer const &) = default;
  ThresholdComputer &operator=(ThresholdComputer const &) = delete;
  ThresholdComputer(ThresholdComputer &&) = default;
  ThresholdComputer &operator=(ThresholdComputer &&) = delete;

  /**
   * @brief Gives the level-dependent threshold for the multiresolution
   * analysis.
   * @param level The level on which the threshold is computed.
   * @return Computed threshold.
   */
  double ThresholdOnLevel(unsigned int const level) const {
    return reference_epsilon_ *
           std::pow(2.0, -1.0 * D_ * (double(maximum_level_) - double(level)));
  }
};

using Thresholder = ThresholdComputer<CC::DIM(), CC::STDO()>;

#endif // THRESHOLD_COMPUTER_H
