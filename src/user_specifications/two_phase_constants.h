/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#ifndef TWO_PHASE_CONSTANTS_H
#define TWO_PHASE_CONSTANTS_H

#include "enums/geometry_settings.h"

namespace GeneralTwoPhaseSettings {
/**
 * Indicates whether information about the convergence of the iterative methods is written to the log-file
 * and terminal or not.
 */
constexpr bool LogConvergenceInformation = false;
/**
 * Indicates whether the number of multiphase nodes is to be written in the log-file and terminal or not.
 */
constexpr bool LogMultiPhaseNodeCount = true;

/**
 * Indicates whether the number of leaves holding a levelset buffer is to be written in the log-file and terminal or not.
 */
constexpr bool LogLevelsetLeafCount = true;
}

namespace GeometryCalculationSettings {
/**
 * Indicates which criteria is used to determine cut cells
 */
constexpr CutCellCriteria CutCellCriteria = CutCellCriteria::ValueBased;

/**
 * Indicates whether a derivative or reconstruction stencil should be used inside the geometry calculator.
 */
constexpr GeometryStencilType GeometryStencilType = GeometryStencilType::Reconstruction;
}

namespace ReinitializationConstants {
/**
 * The pseudo-timestep size to iteratively solve the reinitialization equation.
 * NOTE: It is given as a portion of the cell_size! E.g.: 0.5 results in a pseudo-timestep size of 0.5 * cell_size.
 */
constexpr double Dtau = 0.125;

/**
 * The maximum number of iterations used to solve the reinitialization equation.
 */
constexpr unsigned int MaximumNumberOfIterations = 30; // 30
/**
 * The indicator whether a convergence criteria is applied (true) or a fixed number of iterations (false) is used.
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
}

namespace ExtensionConstants {
/**
 * The pseudo-timestep size to iteratively solve the extension equation.
 * NOTE: It is given as a portion of the cell_size! E.g.: 0.5 results in a pseudo-timestep size of 0.5 * cell_size.
 */
constexpr double Dtau = 0.3;

/**
 * The maximum number of iterations used to solve the extension equation.
 */
constexpr unsigned int MaximumNumberOfIterations = 35;

/**
 * The indicator whether a convergence criteria is applied (true) or a fixed number of iterations (false) is used.
 */
constexpr bool TrackConvergence = true;

/**
 * The maximum residuum allowed if a convergence criteria is applied.
 */
constexpr double MaximumResiduum = 1.0e-3;
}

namespace InterfaceStateTreatmentConstants {
/**
 * The indicator whether we extend the interface quantities or not.
 */
constexpr bool ExtendInterfaceQuantities = true;
}

namespace InterfaceStateExtensionConstants {
/**
 * The pseudo-timestep size to iteratively solve the extension equation.
 * NOTE: It is given as a portion of the cell_size! E.g.: 0.5 results in a pseudo-timestep size of 0.5 * cell_size.
 */
constexpr double Dtau = 0.3;

/**
 * The maximum number of iterations used to solve the extension equation.
 */
constexpr unsigned int MaximumNumberOfIterations = 35;

/**
 * The indicator whether a convergence criteria is applied (true) or a fixed number of iterations (false) is used.
 */
constexpr bool TrackConvergence = true;

/**
 * The maximum residuum allowed if a convergence criteria is applied.
 */
constexpr double MaximumResiduum = 1.0e-3;
}

namespace IterativeInterfaceRiemannSolverConstants {
/**
 * The maximum number of iterations used to solve the interface Riemann problem.
 */
constexpr unsigned int MaximumNumberOfIterations = 100;

/**
 * The maximum residuum allowed for the interface Riemann problem
 */
constexpr double MaximumResiduum = 1.0e-6;
}

#endif // TWO_PHASE_CONSTANTS_H
