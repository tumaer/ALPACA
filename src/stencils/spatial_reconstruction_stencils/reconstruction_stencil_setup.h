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
#ifndef RECONSTRUCTION_STENCIL_SETUP_H
#define RECONSTRUCTION_STENCIL_SETUP_H

#include "user_specifications/stencil_setup.h"
#include "first_order.h"
#include "weno3.h"
#include "weno5.h"
#include "weno-cu6.h"
#include "weno5_z.h"
#include "teno5.h"
#include "weno7.h"
#include "weno9.h"
#include "weno-ao53.h"
#include "weno5_hm.h"
#include "fourth_order_central.h"


/**
 * @brief A namespace to get a ReconstructionStencil type based on a specified constexpr.
 */
namespace ReconstructionStencilSetup {

   /**
    * @brief Function returning the typedef of a ReconstructionStencil based on a constexpr template.
    * 
    * @tparam ReconstructionStencils The constexpr template parameter to specify the exact ReconstructionStencil type.
    */
   template<ReconstructionStencils>
   struct Concretize;

   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::FirstOrder> {
      typedef FirstOrder type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::WENO3> {
      typedef WENO3 type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::FourthOrderCentral> {
      typedef FourthOrderCentral type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::WENO5> {
      typedef WENO5 type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::WENO5Z> {
      typedef WENO5Z type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::WENOAO53> {
      typedef WENOAO53 type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::TENO5> {
      typedef TENO5 type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::WENOCU6> {
      typedef WENOCU6 type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::WENO7> {
      typedef WENO7 type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::WENO9> {
      typedef WENO9 type;
   };
   /**
    * @brief See generic implementation.
    */
   template<>
   struct Concretize<ReconstructionStencils::WENO5HM> {
      typedef WENO5HM type;
   };

}

#endif // RECONSTRUCTION_STENCIL_SETUP_H
