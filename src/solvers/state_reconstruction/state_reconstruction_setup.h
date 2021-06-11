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
#ifndef STATE_RECONSTRUCTION_SETUP_H
#define STATE_RECONSTRUCTION_SETUP_H

#include "user_specifications/equation_settings.h"
#include "user_specifications/state_reconstruction_settings.h"
#include "solvers/state_reconstruction/conservative_state_reconstruction.h"
#include "solvers/state_reconstruction/primitive_state_reconstruction.h"
#include "solvers/state_reconstruction/characteristic_state_reconstruction.h"
#include "solvers/state_reconstruction/gamma_primitive_state_reconstruction.h"
#include "solvers/state_reconstruction/gamma_characteristic_state_reconstruction.h"

static_assert( !( active_equations == EquationSet::Isentropic && state_reconstruction_type == StateReconstructionType::Characteristic ), "Characteristic reconstruction not implemented for isentropic equations!" );

/**
 * @brief A namespace to get a StateReconstruction type based on a specified constexpr.
 */
namespace StateReconstructionSetup {

   namespace EulerNavierStokes {
      /**
       * @brief Function returning the Euler or Navier-Stokes equations state reconstruction matching the type in the template argument.
       * @tparam StateReconstructionType Specification of the state reconstruction type.
       */
      template<StateReconstructionType>
      struct Concretize;

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<StateReconstructionType::Conservative> {
         using type = ConservativeStateReconstruction;
      };

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<StateReconstructionType::Primitive> {
         using type = PrimitiveStateReconstruction;
      };
      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<StateReconstructionType::Characteristic> {
         using type = CharacteristicStateReconstruction;
      };
   }// namespace EulerNavierStokes

   namespace Isentropic {
      /**
       * @brief Function returning the Isentropic state reconstruction matching the type in the template argument.
       * @tparam StateReconstructionType Specification of the state reconstruction type.
       */
      template<StateReconstructionType>
      struct Concretize;

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<StateReconstructionType::Conservative> {
         using type = ConservativeStateReconstruction;
      };

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<StateReconstructionType::Primitive> {
         using type = PrimitiveStateReconstruction;
      };
   }// namespace Isentropic

   namespace GammaModel {
      /**
       * @brief Function returning the Gamma Model state reconstruction matching the type in the template argument.
       * @tparam StateReconstructionType Specification of the state reconstruction type.
       */
      template<StateReconstructionType>
      struct Concretize;

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<StateReconstructionType::Conservative> {
         using type = ConservativeStateReconstruction;
      };

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<StateReconstructionType::Primitive> {
         using type = GammaPrimitiveStateReconstruction;
      };

      /**
       * @brief See generic implementation.
       */
      template<>
      struct Concretize<StateReconstructionType::Characteristic> {
         using type = GammaCharacteristicStateReconstruction;
      };
   }// namespace GammaModel

   /**
   * @brief Function returning the state reconstruction type matching the type in the template argument according to the equation set in the second template parameter.
   * @tparam StateReconstructionType Specification of the state reconstruction type.
   * @tparam EquationSet Specification of the equation(s) beeing solved.
   */
   template<StateReconstructionType S, EquationSet>
   struct Dispatch {
      using type = typename EulerNavierStokes::Concretize<S>::type;
   };

   /**
    * @brief See generic implementation.
    */
   template<StateReconstructionType S>
   struct Dispatch<S, EquationSet::Isentropic> {
      using type = typename Isentropic::Concretize<S>::type;
   };

   /**
    * @brief See generic implementation.
    */
   template<StateReconstructionType S>
   struct Dispatch<S, EquationSet::GammaModel> {
      using type = typename GammaModel::Concretize<S>::type;
   };

   /**
    * @brief Function returning the state reconstruction type for the (globally) selected equation and the given reconstruction template argument.
    * @tparam StateReconstructionType Specification of the state reconstruction type.
    */
   template<StateReconstructionType S>
   struct Concretize {
      using type = typename Dispatch<S, active_equations>::type;
   };

}// namespace StateReconstructionSetup

#endif// STATE_RECONSTRUCTION_SETUP_H
