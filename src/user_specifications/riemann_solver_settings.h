//===------------------ riemann_solver_settings.h -------------------------===//
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
#ifndef RIEMANN_SOLVER_SETTINGS_H
#define RIEMANN_SOLVER_SETTINGS_H

#include "enums/flux_splitting.h"
#include "enums/signal_speed.h"

enum class ConvectiveTermSolvers { FluxSplitting, FiniteVolume };
constexpr ConvectiveTermSolvers convective_term_solver =
    ConvectiveTermSolvers::FluxSplitting;

namespace FluxSplittingSettings {
/* FluxSplitting options are:
 * Roe | LocalLaxFriedrichs | GlobalLaxFriedrichs | Roe_M | LocalLaxFriedrichs_M
 * Roe_M and LocalLaxFriedrichs_M according to \cite Fleischmann20
 */
constexpr FluxSplitting flux_splitting_scheme = FluxSplitting::Roe;

/* Phi in \cite Fleischmann20.
 * Limits the speed of sound in the eigenvalue calculation of Roe-M and LLF-M
 */
constexpr double low_mach_number_limit_factor = 5.0;
} // namespace FluxSplittingSettings

namespace FiniteVolumeSettings {
// RIEMANN_SOLVER
enum class RiemannSolvers { Hllc, Hllc_LM, Hll };
constexpr RiemannSolvers riemann_solver = RiemannSolvers::Hllc;

/* Signal Speed choices for HLL-type solvers are:
 * Einfeldt \cite Einfeldt88 | Davis \cite Davis88 | Toro \cite Toro94 |
 * Arithmetic \cite Coralic14
 */
constexpr SignalSpeed signal_speed_selection = SignalSpeed::Einfeldt;
} // namespace FiniteVolumeSettings

#endif // RIEMANN_SOLVER_SETTINGS_H
