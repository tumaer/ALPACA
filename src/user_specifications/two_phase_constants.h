//===--------------------- two_phase_constants.h --------------------------===//
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
#ifndef TWO_PHASE_CONSTANTS_H
#define TWO_PHASE_CONSTANTS_H

#include "enums/geometry_settings.h"

namespace GeneralTwoPhaseSettings {
/**
 * Indicates whether information about the convergence of the iterative methods
 * is written to the log-file and terminal or not.
 */
constexpr bool LogConvergenceInformation = false;
/**
 * Indicates whether the number of multiphase nodes is to be written in the
 * log-file and terminal or not.
 */
constexpr bool LogMultiPhaseNodeCount = true;

/**
 * Indicates whether the number of leaves holding a levelset buffer is to be
 * written in the log-file and terminal or not.
 */
constexpr bool LogLevelsetLeafCount = true;
} // namespace GeneralTwoPhaseSettings

namespace GeometryCalculationSettings {
/**
 * Indicates which criteria is used to determine cut cells
 */
constexpr CutCellCriteria CutCellCriteria = CutCellCriteria::ValueBased;

/**
 * Indicates whether a derivative or reconstruction stencil should be used
 * inside the geometry calculator.
 */
constexpr GeometryStencilType GeometryStencilType =
    GeometryStencilType::Reconstruction;
} // namespace GeometryCalculationSettings

namespace ReinitializationConstants {
/**
 * The pseudo-timestep size to iteratively solve the reinitialization equation.
 * NOTE: It is given as a portion of the cell_size! E.g.: 0.5 results in a
 * pseudo-timestep size of 0.5 * cell_size.
 */
constexpr double Dtau = 0.125;

/**
 * The maximum number of iterations used to solve the reinitialization equation.
 */
constexpr unsigned int MaximumNumberOfIterations = 30; // 30
/**
 * The indicator whether a convergence criteria is applied (true) or a fixed
 * number of iterations (false) is used.
 */
constexpr bool TrackConvergence = true;

/**
 * The maximum residuum allowed if a convergence criteria is applied.
 */
constexpr double MaximumResiduum = 1.0e-3;

/**
 * Decision whether also cut cells are reinitialized.
 */
static constexpr bool ReinitializeCutCells = true;

/**
 * Decision whether only in last RK stage is reinitialized.
 */
static constexpr bool ReinitializeOnlyInLastRkStage = true;

/**
 * Decision whether reinitialization is performed before or after mixing.
 */
static constexpr bool ReinitializeAfterMixing = true;
} // namespace ReinitializationConstants

namespace ExtensionConstants {
/**
 * The pseudo-timestep size to iteratively solve the extension equation.
 * NOTE: It is given as a portion of the cell_size! E.g.: 0.5 results in a
 * pseudo-timestep size of 0.5 * cell_size.
 */
constexpr double Dtau = 0.3;

/**
 * The maximum number of iterations used to solve the extension equation.
 */
constexpr unsigned int MaximumNumberOfIterations = 35;

/**
 * The indicator whether a convergence criteria is applied (true) or a fixed
 * number of iterations (false) is used.
 */
constexpr bool TrackConvergence = true;

/**
 * The maximum residuum allowed if a convergence criteria is applied.
 */
constexpr double MaximumResiduum = 1.0e-3;
} // namespace ExtensionConstants

namespace InterfaceStateTreatmentConstants {
/**
 * The indicator whether we extend the interface quantities or not.
 */
constexpr bool ExtendInterfaceQuantities = true;
} // namespace InterfaceStateTreatmentConstants

namespace InterfaceStateExtensionConstants {
/**
 * The pseudo-timestep size to iteratively solve the extension equation.
 * NOTE: It is given as a portion of the cell_size! E.g.: 0.5 results in a
 * pseudo-timestep size of 0.5 * cell_size.
 */
constexpr double Dtau = 0.3;

/**
 * The maximum number of iterations used to solve the extension equation.
 */
constexpr unsigned int MaximumNumberOfIterations = 35;

/**
 * The indicator whether a convergence criteria is applied (true) or a fixed
 * number of iterations (false) is used.
 */
constexpr bool TrackConvergence = true;

/**
 * The maximum residuum allowed if a convergence criteria is applied.
 */
constexpr double MaximumResiduum = 1.0e-3;
} // namespace InterfaceStateExtensionConstants

namespace IterativeInterfaceRiemannSolverConstants {
/**
 * The maximum number of iterations used to solve the interface Riemann problem.
 */
constexpr unsigned int MaximumNumberOfIterations = 100;

/**
 * The maximum residuum allowed for the interface Riemann problem
 */
constexpr double MaximumResiduum = 1.0e-6;
} // namespace IterativeInterfaceRiemannSolverConstants

#endif // TWO_PHASE_CONSTANTS_H
