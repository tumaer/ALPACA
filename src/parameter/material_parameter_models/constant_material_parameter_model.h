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
#ifndef CONSTANT_MATERIAL_PARAMETER_MODEL_H
#define CONSTANT_MATERIAL_PARAMETER_MODEL_H

#include "parameter/material_parameter_model.h"

/**
 * @brief The ConstantMaterialParameterModel class implements the specific computation for a material parameter model that gives constant values.
 *        It provides the loop structure to compute the parameter for a single given block. The actual implementation of the model param = const. is done
 *        in the derived class.
 * @note This class serves as a compatability class to provide the same behavior for all materials in case one uses a parameter model to avoid additional
 *       if else statements calling if a model is used or not. Therefore, if one material uses a model, all materials at least must define a constant model.
 *
 * @tparam DerivedConstantMaterialParameterModel derived constan model class that provides the model computation.
 */
template<typename DerivedConstantMaterialParameterModel>
class ConstantMaterialParameterModel : public MaterialParameterModel {

   // friend model, which effectively implements the computation of the parameter
   friend DerivedConstantMaterialParameterModel;

   /**
    * @brief Executes the actual parameter calculation on the complete block.
    * @param block Block on which the parameter calculation should be carried out (parameter on block as indirect return).
    */
   void DoUpdateParameter( Block & block, double const ) const override {

      // extract the parameter from the block, which should be computed
      double (&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetParameterBuffer( DerivedConstantMaterialParameterModel::parameter_buffer_type_ );

      // Compute the parameter based on constant values
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
               parameter_buffer[i][j][k] = static_cast< DerivedConstantMaterialParameterModel const& >( *this ).ComputeParameter();
            }
         }
      }
   }

   /**
    * @brief Executes the actual parameter calculation on the on the block for the given material up to the interface.
    * @param block Block on which the parameter calculation should be carried out (parameter on block as indirect return).
    * @param interface_tags Tags describing the interface position.
    * @param material_sign Sign of the material for identification on interface tags.
    */
   void DoUpdateParameter( Block & block,
                           double const,
                           std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()],
                           std::int8_t const material_sign ) const override {

      // extract the parameter from the block, which should be computed
      double (&parameter_buffer)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetParameterBuffer( DerivedConstantMaterialParameterModel::parameter_buffer_type_ );

      // Compute the parameter based on constan values
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {

               parameter_buffer[i][j][k] = interface_tags[i][j][k] * material_sign >= 0 ?
                                           static_cast< DerivedConstantMaterialParameterModel const& >( *this ).ComputeParameter() : 0.0;
            } //k
         } //j
      }//i
   }

   /**
    * @brief Protected constructor to create the material constant model.
    * @note Can only be called from derived classes.
    */
   explicit ConstantMaterialParameterModel() : MaterialParameterModel() {
      /** Empty besides base class constructor call */
   }

public:
   virtual ~ConstantMaterialParameterModel() = default;
   ConstantMaterialParameterModel( ConstantMaterialParameterModel const& ) = delete;
   ConstantMaterialParameterModel& operator=( ConstantMaterialParameterModel const& ) = delete;
   ConstantMaterialParameterModel( ConstantMaterialParameterModel&& ) = delete;
   ConstantMaterialParameterModel& operator=( ConstantMaterialParameterModel&& ) = delete;

};

#endif // CONSTANT_MATERIAL_PARAMETER_MODEL_H