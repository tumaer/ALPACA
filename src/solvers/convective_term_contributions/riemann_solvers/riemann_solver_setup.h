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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#ifndef RIEMANN_SOLVER_SETUP_H
#define RIEMANN_SOLVER_SETUP_H

#include "user_specifications/numerical_setup.h"
#include "user_specifications/equation_settings.h"
#include "solvers/convective_term_contributions/riemann_solvers/hllc_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/hllc_lm_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/isentropic_hllc_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/hll_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/isentropic_hll_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/gamma_hllc_riemann_solver.h"
#include "solvers/convective_term_contributions/riemann_solvers/gamma_hllc_lm_riemann_solver.h"
#include "user_specifications/riemann_solver_settings.h"

/**
 * @brief A namespace to get a RiemannSolver type based on a specified constexpr.
 */
namespace RiemannSolverSetup {

   namespace Isentropic {
      /**
       * @brief Function returning the Isentropic Riemann solver matching the type in the template argument.
       * @tparam RiemannSolvers Specification of the RiemannSolver type.
       */
      template<FiniteVolumeSettings::RiemannSolvers>
      struct Concretize;

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc> {
         using type = IsentropicHllcRiemannSolver;
      };
      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hll> {
         using type = IsentropicHllRiemannSolver;
      };
   }// namespace Isentropic

   namespace EulerNavierStokes {
      /**
       * @brief Function returning the Euler or Navier-Stokes equations Riemann solver matching the type in the template argument.
       * @tparam RiemannSolvers Specification of the RiemannSolver type.
       */
      template<FiniteVolumeSettings::RiemannSolvers>
      struct Concretize;

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc> {
         using type = HllcRiemannSolver;
      };
      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc_LM> {
         using type = HllcLMRiemannSolver;
      };
      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hll> {
         using type = HllRiemannSolver;
      };
   }// namespace EulerNavierStokes

   namespace GammaModel {
      /**
       * @brief Function returning the Gamma Model Riemann solver matching the type in the template argument.
       * @tparam RiemannSolvers Specification of the RiemannSolver type.
       */
      template<FiniteVolumeSettings::RiemannSolvers>
      struct Concretize;

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc> {
         using type = GammaHllcRiemannSolver;
      };
      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<FiniteVolumeSettings::RiemannSolvers::Hllc_LM> {
         using type = GammaHllcLMRiemannSolver;
      };
   }// namespace GammaModel

   /**
   * @brief Function returning the Riemann solver matching the type in the template argument accroding to the equation set in the second template parameter.
   * @tparam RiemannSolvers Specification of the RiemannSolver type.
   * @tparam EquationSet Specification of the Equation(s) beeing solved.
   */
   template<FiniteVolumeSettings::RiemannSolvers R, EquationSet>
   struct Dispatch {
      using type = typename EulerNavierStokes::Concretize<R>::type;
   };

   /**
    * @brief See generic implementation.
    */
   template<FiniteVolumeSettings::RiemannSolvers R>
   struct Dispatch<R, EquationSet::Isentropic> {
      using type = typename Isentropic::Concretize<R>::type;
   };

   /**
    * @brief See generic implementation.
    */
   template<FiniteVolumeSettings::RiemannSolvers R>
   struct Dispatch<R, EquationSet::GammaModel> {
      using type = typename GammaModel::Concretize<R>::type;
   };

   /**
    * @brief Function returning the Riemann solver for the (globally) selected Equation and the given Solver template argument.
    * @tparam RiemannSolvers Specification of the RiemannSolver type.
    */
   template<FiniteVolumeSettings::RiemannSolvers R>
   struct Concretize {
      using type = typename Dispatch<R, active_equations>::type;
   };

}// namespace RiemannSolverSetup

#endif// RIEMANN_SOLVER_SETUP_H
