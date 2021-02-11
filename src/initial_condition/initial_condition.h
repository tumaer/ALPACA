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
#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

#include <memory>

#include "topology/node_id_type.h"
#include "materials/material_definitions.h"
#include "prime_state_initializer.h"
#include "levelset_initializer.h"

/**
 * @brief The InitialCondition class is used to set the state of all cells according to the user input at the beginning of the simulation. It serves as a
 *        simple handler deligating the work to the appropriate members.
´ */
class InitialCondition {

   std::unique_ptr<PrimeStateInitializer const> prime_state_initializer_;
   // for an arbitrary number of levelsets, this needs to be replaced by a vector of unique pointers. Cannot be const due to non-const member functions.
   std::unique_ptr<LevelsetInitializer> levelset_initializer_;

public:
   InitialCondition() = delete;
   explicit InitialCondition( std::unique_ptr<PrimeStateInitializer const> prime_state_initializer,
                              std::unique_ptr<LevelsetInitializer> levelset_initializer ) : prime_state_initializer_( std::move( prime_state_initializer ) ),
                                                                                            levelset_initializer_( std::move( levelset_initializer ) ) {
      /** Empty besides initializer list */
   }
   ~InitialCondition()                         = default;
   InitialCondition( InitialCondition const& ) = delete;
   InitialCondition& operator=( InitialCondition const& ) = delete;
   InitialCondition( InitialCondition&& )                 = delete;
   InitialCondition& operator=( InitialCondition&& ) = delete;

   // Proxy functions to be called from outside
   void GetInitialPrimeStates( nid_t const node_id, MaterialName const material, double ( &prime_state_buffer )[MF::ANOP()][CC::ICX()][CC::ICY()][CC::ICZ()] ) const {
      prime_state_initializer_->GetInitialPrimeStates( node_id, material, prime_state_buffer );
   };
   void GetInitialLevelset( nid_t const node_id, double ( &levelset_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] ) {
      levelset_initializer_->GetInitialLevelset( node_id, levelset_buffer );
   }
   std::vector<MaterialName> GetInitialMaterials( nid_t const node_id ) {
      return levelset_initializer_->GetInitialMaterials( node_id );
   }
};

#endif// INITIAL_CONDITION_H
