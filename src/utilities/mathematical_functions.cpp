//===------------------- mathematical_functions.cpp -----------------------===//
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
#include "utilities/mathematical_functions.h"
#include "user_specifications/compile_time_constants.h"
#include <algorithm>
#include <cmath>

/**
 * @brief Calculate the Godunov Hamiltonian ( GH ) as described in \cite
 * Min2010.
 * @param derivatives A reference to the single components ( level-set
 * derivatives ) for which the GH has to be calculated.
 * @param old_levelset_sign The sign of the original level-set value.
 * @return The GH as a double value.
 */
double GodunovHamiltonian(double const (&derivatives)[DTI(CC::DIM())][2],
                          double const old_levelset_sign) {
  std::array<double, DTI(CC::DIM())> godunov_hamiltonian_contributions;
  for (unsigned int d = 0; d < DTI(CC::DIM()); ++d) {
    godunov_hamiltonian_contributions[d] =
        std::max(std::max(0.0, old_levelset_sign * derivatives[d][0]) *
                     std::max(0.0, old_levelset_sign * derivatives[d][0]),
                 std::min(0.0, old_levelset_sign * derivatives[d][1]) *
                     std::min(0.0, old_levelset_sign * derivatives[d][1]));
  }
  return std::sqrt(ConsistencyManagedSum(godunov_hamiltonian_contributions));
}
