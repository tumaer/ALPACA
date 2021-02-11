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
#ifndef MATERIAL_PARAMETER_MODEL_H
#define MATERIAL_PARAMETER_MODEL_H

#include "block_definitions/block.h"
#include "materials/material_definitions.h"

/**
 * @brief The MaterialPropertyModel class defines an interface for different discretizations of a model to compute material and material pairing parameters
 *        depending on other materials or physical parameter. It works a proxy to define specific models for different material properties.
 */
class MaterialParameterModel {

protected:
   // required functions needed in derived classes
   virtual void DoUpdateParameter( Block& block, double const cell_size ) const = 0;
   virtual void DoUpdateParameter( Block& block,
                                   double const cell_size,
                                   std::int8_t const ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                   std::int8_t const material_sign ) const      = 0;

   // protected default constructor (can only be called from derived classes)
   explicit MaterialParameterModel() = default;

public:
   virtual ~MaterialParameterModel()                       = default;
   MaterialParameterModel( MaterialParameterModel const& ) = delete;
   MaterialParameterModel& operator=( MaterialParameterModel const& ) = delete;
   MaterialParameterModel( MaterialParameterModel&& )                 = delete;
   MaterialParameterModel& operator=( MaterialParameterModel&& ) = delete;

   /**
    * @brief Computes the desired parameter based on the given parameter for a complete block.
    * @param block Block on which parameters are calculated (indirect return).
    * @param cell_size Size of the cell of given block.
    */
   void UpdateParameter( Block& block, double const cell_size ) const {
      DoUpdateParameter( block, cell_size );
   }

   /**
    * @brief Computes the desired parameter based on the given parameter of an interface-containing block.
    * @param mat_block Pair of MaterialName and corresponding block on which parameters are computed (indirect return).
    * @param cell_size Size of the cell of given block.
    * @param interface_tags Interface tags to specify location of the interface on the given block.
    * @param material_sign Sign of the material for identification on interface tags.
    */
   void UpdateParameter( Block& block,
                         double const cell_size,
                         std::int8_t const ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()],
                         std::int8_t const material_sign ) const {
      DoUpdateParameter( block, cell_size, interface_tags, material_sign );
   }
};

#endif//MATERIAL_PARAMETER_MODEL_H
