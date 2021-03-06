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
#ifndef INITIALIZATION_MATERIAL_H
#define INITIALIZATION_MATERIAL_H

#include <unordered_map>
#include <string>
#include <memory>

#include "input_output/input_reader/material_reader/material_reader.h"
#include "materials/material.h"

/**
 * @brief Defines all instantiation functions required for one single material.
 */
namespace Instantiation {

   // Instantiate function for the equation of state
   std::unique_ptr<EquationOfState const> InstantiateEquationOfState( EquationOfStateName const eos_name,
                                                                      std::unordered_map<std::string, double> const& eos_data,
                                                                      UnitHandler const& unit_handler );

   // Instantiate function for the shear viscosity model
   std::unique_ptr<MaterialParameterModel const> InstantiateShearViscosityModel( MaterialPropertyModelName const model_name,
                                                                                 std::unordered_map<std::string, double> const& model_data,
                                                                                 UnitHandler const& unit_handler );

   // Instantiate function for the thermal conductivity model
   std::unique_ptr<MaterialParameterModel const> InstantiateThermalConductivityModel( MaterialPropertyModelName const model_name,
                                                                                      std::unordered_map<std::string, double> const& model_data,
                                                                                      UnitHandler const& unit_handler );

   // instantiate function for the material type
   MaterialType InstantiateMaterialType( unsigned int const material_index, MaterialReader const& material_reader );

   // instantiate function for the complete material
   std::tuple<MaterialType, Material> InstantiateMaterial( unsigned int const material_index,
                                                           MaterialReader const& material_reader,
                                                           UnitHandler const& unit_handler );

}// namespace Instantiation

#endif// INITIALIZATION_MATERIAL_H
