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
#ifndef LEVELSET_REINITIALIZER_H
#define LEVELSET_REINITIALIZER_H


#include "halo_manager.h"
#include "user_specifications/numerical_setup.h"
#include "levelset/geometry/geometry_calculator.h"
#include "materials/material_manager.h"
#include "levelset/geometry/geometry_calculator_marching_cubes.h"

/**
 * @brief The class LevelsetReinitializer ensures the (signed-)distance property of a level-set field.
 * @tparam DerivedLevelsetReinitializer Typename as template parameter due to CRTP.
 */
template<typename DerivedLevelsetReinitializer>
class LevelsetReinitializer {

   friend DerivedLevelsetReinitializer;

protected:
   HaloManager & halo_manager_; // TODO-19 NH Think about making it const (rats tail)
   LogWriter& logger_;

   /**
    * @brief The default constructor for a LevelsetReinitializer object.
    * @param halo_manager Instance to a HaloManager which provides MPI-related methods.
    */
   explicit LevelsetReinitializer( HaloManager & halo_manager ) :
       halo_manager_( halo_manager ),
       logger_( LogWriter::Instance() )
   {
      // Empty Constructor, besides initializer list.
   }

public:
   LevelsetReinitializer() = delete;
   ~LevelsetReinitializer() = default;
   LevelsetReinitializer( LevelsetReinitializer const& ) = delete;
   LevelsetReinitializer& operator=( LevelsetReinitializer const& ) = delete;
   LevelsetReinitializer( LevelsetReinitializer&& ) = delete;
   LevelsetReinitializer& operator=( LevelsetReinitializer&& ) = delete;

   /**
    * @brief The method to reinitialize a level-set field.
    * @param nodes The nodes for which the level-set field is reinitialized.
    * @param stage The current stage of the Runge-Kutta method.
    */
   void Reinitialize(std::vector<std::reference_wrapper<Node>> const& nodes, bool const is_last_stage) const {
      static_cast<DerivedLevelsetReinitializer const&>(*this).ReinitializeImplementation(nodes, is_last_stage);
   }
};


#endif //LEVELSET_REINITIALIZER_H
