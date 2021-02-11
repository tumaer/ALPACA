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
#include <catch2/catch.hpp>
#include <limits>

#include <limits>

#include "materials/equations_of_state/isentropic.h"
#include "materials/equations_of_state/stiffened_gas.h"
#include "materials/equations_of_state/stiffened_gas_safe.h"
#include "materials/equations_of_state/stiffened_gas_complete_safe.h"
#include "materials/equations_of_state/waterlike_fluid.h"
#include "materials/equations_of_state/noble_abel_stiffened_gas.h"

namespace TestEos {
   /**
    * @brief Tests the getter functions of an equation of state.
    * @param equation_of_state The equation of state under test.
    * @param expected_gamma The expected gamma.
    * @param expected_B The expected B value.
    * @param expected_gruneisen The expected value for the Grunesien parameter.
    * @param gruneisen_is_density_dependent Indicator, whether a density needs to be supplied to the Gruneisen getter.
    * @param comparison_epsilon The prescision needed to be fullfiled in the comparions.
    * @note If Gruneisen is density dependent a density of 42.77 is used.
    */
   void GetterFunctions( std::unique_ptr<EquationOfState const> const& equation_of_state, double const expected_gamma, double const expected_B, double const expected_gruneisen, bool const gruneisen_is_density_dependent = false, double const comparison_epsilon = std::numeric_limits<double>::epsilon() ) {
      REQUIRE( equation_of_state->Gamma() == expected_gamma );
      REQUIRE( equation_of_state->B() == expected_B );
      if( gruneisen_is_density_dependent ) {
         REQUIRE( equation_of_state->Gruneisen( 42.77 ) == Approx( expected_gruneisen ).epsilon( comparison_epsilon ) );
      } else {
         REQUIRE( equation_of_state->Gruneisen() == Approx( expected_gruneisen ).epsilon( comparison_epsilon ) );
      }
   }

   /**
    * @brief Tests the computation of quantities (prime states) from conservatives inputs for the given equation of state.
    * @param equation_of_state The equation of state under test.
    * @param expected_pressure The expected pressure.
    * @param expected_enthalpy The expected enthalpy.
    * @param expected_temperature The expected temperature.
    * @param mass The mass value used in the calculations.
    * @param energy The energy value used in the calculations.
    * @param momentum_x The momentum value used in the calculations.
    * @param momentum_y The momentum value used in the calculations.
    * @param momentum_z The momentum value used in the calculations.
    * @param comparison_epsilon The prescision needed to be fullfiled in the comparions.
    */
   void CalculationsFromConservatives( std::unique_ptr<EquationOfState const> const& equation_of_state, double const expected_pressure, double const expected_enthalpy, double const expected_temperature, double const mass, double const energy, double const momentum_x, double const momentum_y, double const momentum_z, double const comparison_epsilon = std::numeric_limits<double>::epsilon() ) {
      REQUIRE( equation_of_state->Pressure( mass, momentum_x, momentum_y, momentum_z, energy ) == Approx( expected_pressure ).epsilon( comparison_epsilon ) );
      REQUIRE( equation_of_state->Enthalpy( mass, momentum_x, momentum_y, momentum_z, energy ) == Approx( expected_enthalpy ).epsilon( comparison_epsilon ) );
      REQUIRE( equation_of_state->Temperature( mass, momentum_x, momentum_y, momentum_z, energy ) == Approx( expected_temperature ).epsilon( comparison_epsilon ) );
   }

   /**
    * @brief Tests the computation of quantities (conservatives and prime states) from prime state inputs for the given equation of state.
    * @param equation_of_state The equation of state under test.
    * @param expected_energy The expected energy.
    * @param expected_speed_of_sound The expected speed of sound.
    * @param expected_psi The expected psi.
    * @param density The density value used in the calculations.
    * @param pressure The pressure value used in the calculations.
    * @param velocity_x The velocity value used in the calculations.
    * @param velocity_y The velocity value used in the calculations.
    * @param velocity_z The velocity value used in the calculations.
    * @param comparison_epsilon The prescision needed to be fullfiled in the comparions.
    */
   void CalculationsFromPrimes( std::unique_ptr<EquationOfState const> const& equation_of_state, double const expected_energy, double const expected_speed_of_sound, double const expected_psi, double const density, double const pressure, double const velocity_x, double const velocity_y, double const velocity_z, double const comparison_epsilon = std::numeric_limits<double>::epsilon() ) {
      REQUIRE( equation_of_state->Energy( density, velocity_x, velocity_y, velocity_z, pressure ) == Approx( expected_energy ).epsilon( comparison_epsilon ) );
      REQUIRE( equation_of_state->SpeedOfSound( density, pressure ) == Approx( expected_speed_of_sound ).epsilon( comparison_epsilon ) );
      REQUIRE( equation_of_state->Psi( pressure, 1.0 / density ) == Approx( expected_psi ).epsilon( comparison_epsilon ) );
   }

   /**
    * @brief Proviedes some conservative values.
    * @return Mass, energy, and x-, y-, z-momentum.
    */
   constexpr auto ArbitraryConservatives() {
      return std::make_tuple( 1.2, 4.5, 0.36, 2.7, 3.05 );
   }

   /**
    * @brief Proviedes some prime state values.
    * @return Density, pressure, and x-, y-, z-velocity.
    */
   constexpr auto ArbitraryPrimes() {
      return std::make_tuple( 2.4, 3.75, 0.75, 1.9, 4.003 );
   }
}// namespace TestEos

SCENARIO( "Equations of state getters and computation", "[1rank]" ) {
   // Initialize the unit handler class
   UnitHandler const unit_handler( 1.0, 1.0, 1.0, 1.0 );

   GIVEN( "A stiffened gas with a set gamma and background pressure" ) {
      constexpr double gamma                                 = 1.6;
      constexpr double background_pressure                   = 2.0;
      std::unordered_map<std::string, double> const eos_data = { { "gamma", gamma }, { "backgroundPressure", background_pressure } };
      std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGas const>( eos_data, unit_handler ) );

      WHEN( "We check the set members" ) {
         THEN( "The values correspond to the given input" ) {
            TestEos::GetterFunctions( equation_of_state, gamma, background_pressure, 0.6 );
         }
      }

      WHEN( "We set some arbitrary conservative and prime state values" ) {
         // As of now C++ does not allow constexpr decompositions, but it might change in the future (see follow-up of P0488R0)
         auto const [mass, energy, momentum_x, momentum_y, momentum_z]      = TestEos::ArbitraryConservatives();
         auto const [density, pressure, velocity_x, velocity_y, velocity_z] = TestEos::ArbitraryPrimes();
         THEN( "The quantities computed from conservatives match the expectation" ) {
            TestEos::CalculationsFromConservatives( equation_of_state, -4.6805249999999994, -0.1504374999999995, -1, mass, energy, momentum_x, momentum_y, momentum_z );
         }
         THEN( "The quantities computed from primes" ) {
            TestEos::CalculationsFromPrimes( equation_of_state, 35.819144133333332, 1.957890020745122, 2.8958333333333335, density, pressure, velocity_x, velocity_y, velocity_z );
         }
      }
   }

   GIVEN( "A stiffened gas complete safe with a set gamma and background pressure, thermal energy factor and specific gas constant" ) {
      constexpr double gamma                                 = 4.4;
      constexpr double background_pressure                   = 6e8;
      constexpr double thermal_energy_factor                 = 6.5;
      constexpr double energy_translation_factor             = 0.76;
      constexpr double specific_gas_constant                 = 1.4;
      std::unordered_map<std::string, double> const eos_data = { { "gamma", gamma }, { "energyTranslationFactor", energy_translation_factor }, { "backgroundPressure", background_pressure }, { "thermalEnergyFactor", thermal_energy_factor }, { "specificGasConstant", specific_gas_constant } };
      std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGasCompleteSafe const>( eos_data, unit_handler ) );

      WHEN( "We check the set members" ) {
         THEN( "The values correspond to the given input" ) {
            TestEos::GetterFunctions( equation_of_state, gamma, background_pressure, 3.4 );
         }
      }

      WHEN( "We set some arbitrary conservative and prime state values" ) {
         // As of now C++ does not allow constexpr decompositions, but it might change in the future (see follow-up of P0488R0)
         auto const [mass, energy, momentum_x, momentum_y, momentum_z]      = TestEos::ArbitraryConservatives();
         auto const [density, pressure, velocity_x, velocity_y, velocity_z] = TestEos::ArbitraryPrimes();
         THEN( "The quantities computed from conservatives match the expectation" ) {
            TestEos::CalculationsFromConservatives( equation_of_state, -600000000.0, -499999996.25, 29.341375851961814, mass, energy, momentum_x, momentum_y, momentum_z );
         }
         THEN( "The quantities computed from primes" ) {
            TestEos::CalculationsFromPrimes( equation_of_state, 776470615.39804602, 33166.248007198526, 1100000001.5625, density, pressure, velocity_x, velocity_y, velocity_z );
         }
      }
   }

   GIVEN( "A safe stiffened gas with a set gamma and background pressure" ) {
      constexpr double gamma                                 = 1.2;
      constexpr double background_pressure                   = 0.2;
      std::unordered_map<std::string, double> const eos_data = { { "gamma", gamma }, { "backgroundPressure", background_pressure } };
      std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<StiffenedGasSafe const>( eos_data, unit_handler ) );

      WHEN( "We check the set members" ) {
         THEN( "The values correspond to the given input" ) {
            TestEos::GetterFunctions( equation_of_state, gamma, background_pressure, 0.2 );
         }
      }

      WHEN( "We set some arbitrary conservative and prime state values" ) {
         // As of now C++ does not allow constexpr decompositions, but it might change in the future (see follow-up of P0488R0)
         auto const [mass, energy, momentum_x, momentum_y, momentum_z]      = TestEos::ArbitraryConservatives();
         auto const [density, pressure, velocity_x, velocity_y, velocity_z] = TestEos::ArbitraryPrimes();
         THEN( "The quantities computed from conservatives match the expectation" ) {
            TestEos::CalculationsFromConservatives( equation_of_state, -0.19999999999999979, 3.5833333333333335, -1.0, mass, energy, momentum_x, momentum_y, momentum_z );
         }
         THEN( "The quantities computed from primes" ) {
            TestEos::CalculationsFromPrimes( equation_of_state, 44.185810800000006, 1.4053469322555197, 1.6625000000000001, density, pressure, velocity_x, velocity_y, velocity_z );
         }
      }
   }

   GIVEN( "A waterlike fluid with set gamma, A, B, and rho0" ) {
      constexpr double gamma                                 = 6.8;
      constexpr double A                                     = 2.5;
      constexpr double B                                     = 6.8;
      constexpr double rho0                                  = 4.4;
      std::unordered_map<std::string, double> const eos_data = { { "gamma", gamma }, { "A", A }, { "B", B }, { "rho0", rho0 } };
      std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<WaterlikeFluid const>( eos_data, unit_handler ) );

      WHEN( "We check the set members" ) {
         THEN( "The values correspond to the given input" ) {
            TestEos::GetterFunctions( equation_of_state, -1.0, -1, 0.0 );
         }
      }

      WHEN( "We set some arbitrary conservative and prime state values" ) {
         // As of now C++ does not allow constexpr decompositions, but it might change in the future (see follow-up of P0488R0)
         auto const [mass, energy, momentum_x, momentum_y, momentum_z]      = TestEos::ArbitraryConservatives();
         auto const [density, pressure, velocity_x, velocity_y, velocity_z] = TestEos::ArbitraryPrimes();
         THEN( "The quantities computed from conservatives match the expectation" ) {
            TestEos::CalculationsFromConservatives( equation_of_state, -4.2990103920267781, 0.0, -1.0, mass, energy, momentum_x, momentum_y, momentum_z );
         }
         THEN( "The quantities computed from primes" ) {
            TestEos::CalculationsFromPrimes( equation_of_state, 27.147879765517239, 4.7758070871145257, 22.808333333333334, density, pressure, velocity_x, velocity_y, velocity_z );
         }
      }
   }

   GIVEN( "A Noble-Abel stiffened gas with set gamma, covolume, pressure_constant, energy_constant, entropy_constant and specific_heat_capacity" ) {
      constexpr double gamma                                 = 6.8;
      constexpr double covolume                              = 0.04;
      constexpr double pressure_constant                     = 4.3;
      constexpr double energy_constant                       = 2.0;
      constexpr double entropy_constant                      = 2.0;
      constexpr double specific_heat_capacity                = 4.0;
      std::unordered_map<std::string, double> const eos_data = { { "gamma", gamma }, { "covolume", covolume }, { "pressureConstant", pressure_constant }, { "energyConstant", energy_constant }, { "entropyConstant", entropy_constant }, { "specificHeatCapacity", specific_heat_capacity } };
      std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<NobleAbelStiffenedGas const>( eos_data, unit_handler ) );

      WHEN( "We check the set members" ) {
         THEN( "The values correspond to the given input" ) {
            TestEos::GetterFunctions( equation_of_state, gamma, 4.3, -8.1598199212155311, true );
         }
      }

      WHEN( "We set some arbitrary conservative and prime state values" ) {
         // As of now C++ does not allow constexpr decompositions, but it might change in the future (see follow-up of P0488R0)
         auto const [mass, energy, momentum_x, momentum_y, momentum_z]      = TestEos::ArbitraryConservatives();
         auto const [density, pressure, velocity_x, velocity_y, velocity_z] = TestEos::ArbitraryPrimes();
         THEN( "The quantities computed from conservatives match the expectation" ) {
            TestEos::CalculationsFromConservatives( equation_of_state, -58.89519082633052, -45.329325688608769, -1.8669045138888885, mass, energy, momentum_x, momentum_y, momentum_z );
         }
         THEN( "The quantities computed from primes" ) {
            TestEos::CalculationsFromPrimes( equation_of_state, 34.177700455172413, 5.0229928555731238, 15.205567846607668, density, pressure, velocity_x, velocity_y, velocity_z );
         }
      }
   }

   GIVEN( "An isentropic equation of state with set gamma and A" ) {
      constexpr double gamma                                 = 3.2;
      constexpr double A                                     = 0.74586;
      std::unordered_map<std::string, double> const eos_data = { { "gamma", gamma }, { "A", A } };
      std::unique_ptr<EquationOfState const> equation_of_state( std::make_unique<Isentropic const>( eos_data, unit_handler ) );

      WHEN( "We check the set members" ) {
         THEN( "The values correspond to the given input" ) {
            TestEos::GetterFunctions( equation_of_state, gamma, -1.0, -1.0 );
         }
      }

      WHEN( "We set some arbitrary conservative and prime state values" ) {
         // As of now C++ does not allow constexpr decompositions, but it might change in the future (see follow-up of P0488R0)
         auto const [mass, energy, momentum_x, momentum_y, momentum_z]      = TestEos::ArbitraryConservatives();
         auto const [density, pressure, velocity_x, velocity_y, velocity_z] = TestEos::ArbitraryPrimes();
         THEN( "The quantities computed from conservatives match the expectation" ) {
            TestEos::CalculationsFromConservatives( equation_of_state, 1.3367103297833645, -1.0, -1.0, mass, energy, momentum_x, momentum_y, momentum_z );
         }
         THEN( "The quantities computed from primes" ) {
            TestEos::CalculationsFromPrimes( equation_of_state, -1.0, 2.2360679774997898, -1.0, density, pressure, velocity_x, velocity_y, velocity_z );
         }
      }
   }
}
