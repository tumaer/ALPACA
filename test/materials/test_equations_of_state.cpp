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
#include <catch.hpp>

#include "materials/equations_of_state/stiffened_gas.h"
#include "materials/equations_of_state/stiffened_gas_safe.h"
#include "materials/equations_of_state/stiffened_gas_complete_safe.h"
#include "materials/equations_of_state/waterlike_fluid.h"
#include "materials/equations_of_state/noble_abel_stiffened_gas.h"

SCENARIO( "Equations of state are constructed and their properties are queried", "[1rank]" ) {

   // Initialize the unit handler class
   UnitHandler const unit_handler( 1.0, 1.0, 1.0, 1.0 );

   GIVEN( "A stiffened gas" ) {

      WHEN( "All parameters are set to zero" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 0.0 }, { "backgroundPressure", 0.0 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGas const>( eos_data, unit_handler ) );

         THEN( "The temperature is" ) {
            REQUIRE( equation_of_state->GetTemperature( 0.0, 0.0, 0.0, 0.0, 0.0 ) == -1 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen() == -1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( equation_of_state->GetPsi( 0.0, 0.0 ) == 0 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == 0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == 0 );
         }
      }

      WHEN( "The parameters are reasonable values" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 6.1 }, { "backgroundPressure", 4.3 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGas const>( eos_data, unit_handler ) );

         THEN( "The pressure is" ) {
            REQUIRE( equation_of_state->GetPressure( 1.0, 2.0, 3.0, 4.0, 5.0 ) == Approx( -74.68 ) );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( equation_of_state->GetEnthalpy( 6.0, 7.0, 8.0, 9.0, 10.0 ) == Approx( -7.9466666667 ) );
         }

         THEN( "The energy is" ) {
            REQUIRE( equation_of_state->GetEnergy( 11.0, 12.0, 13.0, 14.0, 15.0 ) == Approx( 31.2206773619 ) );
         }

         THEN( "The temperature is" ) {
            REQUIRE( equation_of_state->GetTemperature( 16.0, 17.0, 18.0, 19.0, 20.0 ) == -1 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen() == 5.1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( equation_of_state->GetPsi( 21.0, 22.0 ) == 1039.06 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == 6.1 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == 4.3 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( equation_of_state->GetSpeedOfSound( 23.0, 24.0 ) == Approx( 2.7396445342 ) );
         }
      }
   }


   GIVEN( "A stiffened gas complete safe" ) {

      WHEN( "All parameters are set to zero" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 0.0 }, { "energyTranslationFactor", 0.0 }, { "backgroundPressure", 0.0 }, { "thermalEnergyFactor", 0.0 }, { "specificGasConstant", 0.0 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGasCompleteSafe const>( eos_data, unit_handler ) );

         THEN( "The pressure is" ) {
            REQUIRE( equation_of_state->GetPressure( 0.0, 0.0, 0.0, 0.0, 0.0 ) == Approx( 0.0 ).margin( 0.000001 ) );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( equation_of_state->GetEnthalpy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 1.0 );
         }

         THEN( "The energy is" ) {
            REQUIRE( equation_of_state->GetEnergy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 0.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen() == -1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( equation_of_state->GetPsi( 0.0, 0.0 ) == 0 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == 0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == 0 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( equation_of_state->GetSpeedOfSound( 0.0, 0.0 ) == Approx( 0.0000000149 ).epsilon( 0.1 ) );
         }
      }

      WHEN( "The parameters contain reasonable values" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 6.1 }, { "energyTranslationFactor", 1.0 }, { "backgroundPressure", 4.3 }, { "thermalEnergyFactor", 2.0 }, { "specificGasConstant", 3.0 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGasCompleteSafe const>( eos_data, unit_handler ) );

         THEN( "The pressure is" ) {
            REQUIRE( equation_of_state->GetPressure( 1.0, 2.0, 3.0, 4.0, 5.0 ) == -4.3 );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( equation_of_state->GetEnthalpy( 6.0, 7.0, 8.0, 9.0, 10.0 ) == Approx( 0.95 ) );
         }

         THEN( "The energy is" ) {
            REQUIRE( equation_of_state->GetEnergy( 11.0, 12.0, 13.0, 14.0, 15.0 ) == Approx( 42.2206773619 ) );
         }

         THEN( "The temperature is" ) {
            REQUIRE( equation_of_state->GetTemperature( 16.0, 17.0, 18.0, 19.0, 20.0 ) == Approx( 4704254.7119584298 ) );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen() == 5.1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( equation_of_state->GetPsi( 21.0, 22.0 ) == 1039.06 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == 6.1 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == 4.3 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( equation_of_state->GetSpeedOfSound( 23.0, 24.0 ) == Approx( 2.7396445342 ) );
         }
      }
   }

   GIVEN( "A stiffened gas safe" ) {

      WHEN( "All parameters are set to zero" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 0.0 }, { "backgroundPressure", 0.0 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGasSafe const>( eos_data, unit_handler ) );

         THEN( "The pressure is" ) {
            REQUIRE( equation_of_state->GetPressure( 0.0, 0.0, 0.0, 0.0, 0.0 ) == Approx( 0.0 ).margin( 0.000001 ) );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( equation_of_state->GetEnthalpy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 1.0 );
         }

         THEN( "The energy is" ) {
            REQUIRE( equation_of_state->GetEnergy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 0.0 );
         }

         THEN( "The temperature is" ) {
            REQUIRE( equation_of_state->GetTemperature( 0.0, 0.0, 0.0, 0.0, 0.0 ) == -1.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen() == -1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( equation_of_state->GetPsi( 0.0, 0.0 ) == 0 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == 0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == 0 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( equation_of_state->GetSpeedOfSound( 0.0, 0.0 ) == Approx( 0.0000000149 ).epsilon( 0.1 ) );
         }
      }

      WHEN( "The parameters contain reasonable values" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 6.1 }, { "backgroundPressure", 4.3 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGasSafe const>( eos_data, unit_handler ) );

         THEN( "The pressure is" ) {
            REQUIRE( equation_of_state->GetPressure( 1.0, 2.0, 3.0, 4.0, 5.0 ) == -4.3 );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( equation_of_state->GetEnthalpy( 6.0, 7.0, 8.0, 9.0, 10.0 ) == Approx( 0.95 ) );
         }

         THEN( "The energy is" ) {
            REQUIRE( equation_of_state->GetEnergy( 11.0, 12.0, 13.0, 14.0, 15.0 ) == Approx( 31.2206773619 ) );
         }

         THEN( "The temperature is" ) {
            REQUIRE( equation_of_state->GetTemperature( 16.0, 17.0, 18.0, 19.0, 20.0 ) == -1.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen() == 5.1 );
         }

         THEN( "The psi is" ) {
            REQUIRE( equation_of_state->GetPsi( 21.0, 22.0 ) == 1039.06 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == 6.1 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == 4.3 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( equation_of_state->GetSpeedOfSound( 23.0, 24.0 ) == Approx( 2.7396445342 ) );
         }
      }
   }


   GIVEN( "A waterlike fluid" ) {

      WHEN( "All parameters are set to zero" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 0.0 }, { "A", 0.0 }, { "B", 0.0 }, { "rho0", 0.0 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<WaterlikeFluid const>( eos_data, unit_handler ) );

         THEN( "The enthalpy is" ) {
            REQUIRE( equation_of_state->GetEnthalpy( 0.0, 0.0, 0.0, 0.0, 0.0 ) == 0.0 );
         }

         THEN( "The temperature is" ) {
            REQUIRE( equation_of_state->GetTemperature( 0.0, 0.0, 0.0, 0.0, 0.0 ) == -1.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen() == 0.0 );
         }

         THEN( "The psi is" ) {
            REQUIRE( equation_of_state->GetPsi( 0.0, 0.0 ) == 0 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == -1.0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == -1.0 );
         }
      }

      WHEN( "The parameters contain reasonable values" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 6.1 }, { "A", 1.0 }, { "B", 4.3 }, { "rho0", 1.0 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<WaterlikeFluid const>( eos_data, unit_handler ) );

         THEN( "The pressure is" ) {
            REQUIRE( equation_of_state->GetPressure( 1.0, 2.0, 3.0, 4.0, 5.0 ) == 1.0 );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( equation_of_state->GetEnthalpy( 6.0, 7.0, 8.0, 9.0, 10.0 ) == 0.0 );
         }

         THEN( "The energy is" ) {
            REQUIRE( equation_of_state->GetEnergy( 11.0, 12.0, 13.0, 14.0, 15.0 ) == Approx( 22.8481283422 ) );
         }

         THEN( "The temperature is" ) {
            REQUIRE( equation_of_state->GetTemperature( 16.0, 17.0, 18.0, 19.0, 20.0 ) == -1.0 );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen() == 0.0 );
         }

         THEN( "The psi is" ) {
            REQUIRE( equation_of_state->GetPsi( 21.0, 22.0 ) == 3261.06 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == -1.0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == -1.0 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( equation_of_state->GetSpeedOfSound( 23.0, 24.0 ) == Approx( 2.690805601 ) );
         }
      }
   }

   GIVEN( "A Noble-Abel stiffened gas" ) {

      WHEN( "All parameters are set to zero" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 0.0 }, { "covolume", 0.0 }, { "pressureConstant", 0.0 }, { "energyConstant", 0.0 },
                                                                    { "entropyConstant", 0.0 }, { "specificHeatCapacity", 0.0 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<NobleAbelStiffenedGas const>( eos_data, unit_handler ) );

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen( 0.0 ) == -1.0 );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == 0.0 );
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == 0.0 );
         }
      }

      WHEN( "The parameters contain reasonable values" ) {

         // Define material properties and initialize material and subclasses properly
         std::unordered_map<std::string, double> const eos_data = { { "gamma", 6.1 }, { "covolume", 0.04 }, { "pressureConstant", 4.3 }, { "energyConstant", 2.0 },
                                                                    { "entropyConstant", 3.0 }, { "specificHeatCapacity", 4.0 } };
         std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<NobleAbelStiffenedGas const>( eos_data, unit_handler ) );

         THEN( "The pressure is" ) {
            REQUIRE( equation_of_state->GetPressure( 1.0, 2.0, 3.0, 4.0, 5.0 ) == Approx( -87.32375 ) );
         }

         THEN( "The enthalpy is" ) {
            REQUIRE( equation_of_state->GetEnthalpy( 6.0, 7.0, 8.0, 9.0, 10.0 ) == Approx( -23.02298245 ) );
         }

         THEN( "The energy is" ) {
            REQUIRE( equation_of_state->GetEnergy( 11.0, 12.0, 13.0, 14.0, 15.0 ) == Approx( 49.66357932 ) );
         }

         THEN( "The temperature is" ) {
            REQUIRE( equation_of_state->GetTemperature( 16.0, 17.0, 18.0, 19.0, 20.0 ) == Approx( -0.6872734375 ) );
         }

         THEN( "The gruneisen is" ) {
            REQUIRE( equation_of_state->GetGruneisen( 20.5 ) == Approx( 28.3333333333333333 ) );
         }

         THEN( "The psi is" ) {
            REQUIRE( equation_of_state->GetPsi( 21.0, 22.0 ) == Approx( 1040.952641 ) );
         }

         THEN( "The gamma is" ) {
            REQUIRE( equation_of_state->GetGamma() == 6.1);
         }

         THEN( "The background pressure is" ) {
            REQUIRE( equation_of_state->GetB() == 4.3 );
         }

         THEN( "The speed of sound is" ) {
            REQUIRE( equation_of_state->GetSpeedOfSound( 23.0, 24.0 ) == Approx( 9.686106141 ) );
         }
      }
   }
}