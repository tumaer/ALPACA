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
#ifndef UNIT_HANDLER_H
#define UNIT_HANDLER_H

#include <vector>
#include <cmath>
#include <array>
#include "log_writer.h"
#include "enums/unit_type.h"

 /** 
 * @brief The UnitHandler class takes care of the conversions between unitless and non-unitless representations of quantities.
 *        Within the Simulation Kernel only non-dimensional quantities are used. In the user input and output,
 *        however, values are given in unit repesentation.
 */
class UnitHandler {
   // Input file reference values and SI units
   double const density_reference_; // input file 
   double const velocity_reference_;// input file 
   double const length_reference_; // inputfile + SI unit [m]
   double const temperature_reference_; // input file + SI unit [K]
   double const time_reference_; // SI unit [s]
   double const mass_reference_; // SI unit [kg]
   // derived reference values for each exisitng unit type 
   double const momentum_reference_;
   double const pressure_reference_;
   double const energy_reference_;
   double const viscosity_reference_;
   double const thermal_conductivity_reference_;
   double const surface_tension_coefficient_reference_;

public :
   UnitHandler() = delete;
   explicit UnitHandler( double const density_reference, double const velocity_reference, double const length_reference, double const temperature_reference );
   ~UnitHandler() = default;
   UnitHandler( UnitHandler const& ) = delete;
   UnitHandler& operator=( UnitHandler const& ) = delete;
   UnitHandler( UnitHandler&& ) = delete;
   UnitHandler& operator=( UnitHandler&& ) = delete;

   // Dimensionalization and non-dimensionalization functions for a single unit type
   double NonDimensionalizeValue( double const value, UnitType const unit_type ) const;
   double DimensionalizeValue( double const value, UnitType const unit_type ) const;
   // Dimensionalization and non-dimensionalization for a set of unit types in nominator and denominator
   double NonDimensionalizeValue( double const value, std::vector<UnitType> const& nominator_units, std::vector<UnitType> const& denominator_units ) const;
   double DimensionalizeValue( double const value, std::vector<UnitType> const& nominator_units, std::vector<UnitType> const& denominator_units ) const;
};

#endif // UNIT_HANDLER_H
