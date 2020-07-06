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
#ifndef CONSTANT_SURFACE_TENSION_COEFFICIENT_MODEL_H
#define CONSTANT_SURFACE_TENSION_COEFFICIENT_MODEL_H

#include "unit_handler.h"
#include "parameter/interface_parameter_models/constant_interface_parameter_model.h"

/**
 * @brief Class that provides the actual model implementation f = const. for a Constant material parameter model to compute
 *        the thermal conductivity in the field.
 */
class ConstantSurfaceTensionCoefficientModel : public ConstantInterfaceParameterModel<ConstantSurfaceTensionCoefficientModel> {

   // Defines the friend to be used for the computation indicating the primestates this model depends on
   friend ConstantInterfaceParameterModel;

   // Definition of parameters needed in the base class for the computations
   static constexpr InterfaceParameter parameter_buffer_type_ = InterfaceParameter::SurfaceTensionCoefficient;

   // member variables needed for the model calculation
   double const sigma_constant_;

   // Actual computation of surface tension coefficient
   double ComputeParameter() const;

public:
   ConstantSurfaceTensionCoefficientModel() = delete;
   explicit ConstantSurfaceTensionCoefficientModel( std::unordered_map<std::string, double> const& parameter_map,
                                              UnitHandler const& unit_handler );
   virtual ~ConstantSurfaceTensionCoefficientModel() = default;
   ConstantSurfaceTensionCoefficientModel( const ConstantSurfaceTensionCoefficientModel& ) = delete;
   ConstantSurfaceTensionCoefficientModel& operator=( const ConstantSurfaceTensionCoefficientModel& ) = delete;
   ConstantSurfaceTensionCoefficientModel( ConstantSurfaceTensionCoefficientModel&& ) = delete;
   ConstantSurfaceTensionCoefficientModel& operator=( ConstantSurfaceTensionCoefficientModel&& ) = delete;

   // logging function
   std::string GetLogData( unsigned int const indent, UnitHandler const& unit_handler ) const;
};

#endif // CONSTANT_SURFACE_TENSION_COEFFICIENT_MODEL_H