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

#ifndef TEMPERATURE_INTERFACE_PARAMETER_MODEL_H
#define TEMPERATURE_INTERFACE_PARAMETER_MODEL_H

#include "parameter/interface_parameter_model.h"

/**
 * @brief The TemperatureInterfaceParameterModel class implements the specific computation for a interface parameter model that gives temperature dependent values.
 *        It provides the loop structure to compute the parameter for a single given node. The actual implementation of the model param = f(T) is done
 *        in the derived class.
 *
 * @tparam DerivedTemperatureInterfaceParameterModel derived temperature model class that provides the model computation.
 */
template<typename DerivedTemperatureInterfaceParameterModel>
class TemperatureInterfaceParameterModel : public InterfaceParameterModel {

   // friend model, which effectively implements the computation of the parameter
   friend DerivedTemperatureInterfaceParameterModel;

   /**
    * @brief Executes the actual parameter calculation on the complete block.
    * @param node Node for which the parameter is computed.
    */
   void DoUpdateParameter( Node& node ) const override {

      // extract the temperature from the block
      double const( &T_positive )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( MaterialSignCapsule::PositiveMaterial() ).GetPrimeStateBuffer( PrimeState::Temperature );
      double const( &T_negative )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( MaterialSignCapsule::NegativeMaterial() ).GetPrimeStateBuffer( PrimeState::Temperature );

      // extract the parameter from the block, which should be computed
      double( &parameter_buffer )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetInterfaceParameterBuffer( DerivedTemperatureInterfaceParameterModel::parameter_buffer_type_ );

      // Obtain the interface tags and volume fraction
      std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
      double const( &volume_fraction )[CC::TCX()][CC::TCY()][CC::TCZ()]     = node.GetInterfaceBlock().GetBaseBuffer( Levelset::VolumeFraction );

      // Compute the parameter based on constan values
      for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {

               // Only compute the parameter for cut cells. Average quantities based on the volume fraction
               if( interface_tags[i][j][k] == ITTI( IT::OldCutCell ) ) {
                  double const averaged_T   = volume_fraction[i][j][k] * T_positive[i][j][k] + ( 1.0 - volume_fraction[i][j][k] ) * T_negative[i][j][k];
                  parameter_buffer[i][j][k] = static_cast<DerivedTemperatureInterfaceParameterModel const&>( *this ).ComputeParameter( averaged_T );
               } else {
                  parameter_buffer[i][j][k] = 0.0;
               }
            }//k
         }   //j
      }      //i
   }

   /**
    * @brief Protected constructor to create the interface temperature model.
    * @note Can only be called from derived classes.
    */
   explicit TemperatureInterfaceParameterModel() : InterfaceParameterModel() {
      /** Empty besides base class constructor call */
   }

public:
   virtual ~TemperatureInterfaceParameterModel()                                   = default;
   TemperatureInterfaceParameterModel( TemperatureInterfaceParameterModel const& ) = delete;
   TemperatureInterfaceParameterModel& operator=( TemperatureInterfaceParameterModel const& ) = delete;
   TemperatureInterfaceParameterModel( TemperatureInterfaceParameterModel&& )                 = delete;
   TemperatureInterfaceParameterModel& operator=( TemperatureInterfaceParameterModel&& ) = delete;
};

#endif// TEMPERATURE_INTERFACE_PARAMETER_MODEL_H