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
#ifndef TEMPERATURE_PRESSURE_MATERIAL_PARAMETER_MODEL_H
#define TEMPERATURE_PRESSURE_MATERIAL_PARAMETER_MODEL_H

#include "parameter/material_parameter_model.h"

/**
 * @brief The TemperaturePressureMaterialParameterModel class implements the specific computation for a material parameter model that acts on the prime states
 *        temperature and pressure. It provides the loop structure to compute the parameter for a single given block. The actual implementation of the
 *        model param = f(T,p) is done in the derived class.
 *
 * @tparam DerivedTemperaturePressureMaterialParameterModel derived temperature pressure model class that provides the model computations.
 */
template<typename DerivedTemperaturePressureMaterialParameterModel>
class TemperaturePressureMaterialParameterModel : public MaterialParameterModel {

   // friend model, which effectively implements the computation of the parameter
   friend DerivedTemperaturePressureMaterialParameterModel;

   /**
    * @brief Executes the actual parameter calculation on the complete block.
    * @param block Block on which the parameter calculation should be carried out (parameter on block as indirect return).
    */
   void DoUpdateParameter( Block& block, double const ) const override {

      // extract the temperature from the block
      double const( &T )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetPrimeStateBuffer( PrimeState::Temperature );
      double const( &p )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetPrimeStateBuffer( PrimeState::Pressure );

      // extract the shear viscosity from the block, which should be computed
      double( &parameter_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetParameterBuffer( DerivedTemperaturePressureMaterialParameterModel::parameter_to_calculate_ );

      // Compute the parameter based on the temperature
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {

               parameter_buffer[i][j][k] = static_cast<DerivedTemperaturePressureMaterialParameterModel const&>( *this ).ComputeParameter( T[i][j][k], p[i][j][k] );
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
   void DoUpdateParameter( Block& block,
                           double const,
                           std::int8_t const ( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()],
                           std::int8_t const material_sign ) const override {

      // extract the temperature and pressure from the block
      double const( &T )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetPrimeStateBuffer( PrimeState::Temperature );
      double const( &p )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetPrimeStateBuffer( PrimeState::Pressure );

      // extract the parameter from the block, which should be computed
      double( &parameter_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetParameterBuffer( DerivedTemperaturePressureMaterialParameterModel::parameter_to_calculate_ );

      // Compute the parameter based on the shear rate and model parameter
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {

               // Note: Volume fraction can be used to average temperature and pressure
               parameter_buffer[i][j][k] = interface_tags[i][j][k] * material_sign >= 0 ?
                                                 static_cast<DerivedTemperaturePressureMaterialParameterModel const&>( *this ).ComputeParameter( T[i][j][k], p[i][j][k] ) :
                                                 0.0;

            }//k
         }   //j
      }      //i
   }

   /**
    * @brief Protected constructor to create the material temperature pressure model.
    * @note Can only be called from derived classes.
    */
   explicit TemperaturePressureMaterialParameterModel() : MaterialParameterModel() {
      /** Empty besides base class constructor call */
   }

public:
   virtual ~TemperaturePressureMaterialParameterModel()                                          = default;
   TemperaturePressureMaterialParameterModel( TemperaturePressureMaterialParameterModel const& ) = delete;
   TemperaturePressureMaterialParameterModel& operator=( TemperaturePressureMaterialParameterModel const& ) = delete;
   TemperaturePressureMaterialParameterModel( TemperaturePressureMaterialParameterModel&& )                 = delete;
   TemperaturePressureMaterialParameterModel& operator=( TemperaturePressureMaterialParameterModel&& ) = delete;
};

#endif// TEMPERATURE_PRESSURE_MATERIAL_PARAMETER_MODEL_H
