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
#ifndef STIFFENED_GAS_COMPLETE_SAFE_H
#define STIFFENED_GAS_COMPLETE_SAFE_H

#include <limits>
#include <unordered_map>
#include "unit_handler.h"
#include "materials/equation_of_state.h"

/**
 * @brief The StiffenedGasCompleteSafe class implements a generic stiffened gas, i.e. material parameters must be set via caller input.
 *        The equation of state is described in \cite Menikoff1989. Safe version, i.e. uses epsilons to avoid floating-point errors.
 */
class StiffenedGasCompleteSafe : public EquationOfState {
   // parameters required for the computataion
   double const epsilon_ = std::numeric_limits<double>::epsilon();
   double const gamma_;
   double const energy_translation_factor_;
   double const stiffened_pressure_constant_;
   double const thermal_energy_factor_;
   double const specific_gas_constant_;

   // functions required from the bae class
   double DoGetPressure   ( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const override;
   double DoGetEnthalpy   ( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const override;
   double DoGetEnergy     ( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const pressure ) const override;
   double DoGetTemperature( double const density, double const momentum_x, double const momentum_y, double const momentum_z, double const energy ) const override;
   double DoGetGruneisen() const override;
   double DoGetPsi( double const pressure, double const one_density ) const override;
   double DoGetGamma() const override;
   double DoGetB() const override;
   double DoGetSpeedOfSound( double const density, double const pressure ) const override;

public:
   StiffenedGasCompleteSafe() = delete;
   explicit StiffenedGasCompleteSafe( std::unordered_map<std::string, double> const& dimensional_eos_data, UnitHandler const& unit_handler );
   virtual ~StiffenedGasCompleteSafe() = default;
   StiffenedGasCompleteSafe( StiffenedGasCompleteSafe const& ) = delete;
   StiffenedGasCompleteSafe& operator=( StiffenedGasCompleteSafe const& ) = delete;
   StiffenedGasCompleteSafe( StiffenedGasCompleteSafe&& ) = delete;
   StiffenedGasCompleteSafe operator=( StiffenedGasCompleteSafe&& ) = delete;

   // function for logging
   std::string GetLogData( unsigned int const indent, UnitHandler const& unit_handler ) const;
};

#endif // STIFFENED_GAS_COMPLETE_SAFE_H
