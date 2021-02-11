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
#ifndef THRESHOLD_COMPUTER_H
#define THRESHOLD_COMPUTER_H

#include <cmath>
#include "enums/dimension_definition.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief A class to provide the multiresolution threshold values for the given intput. According to \cite Harten1995 and \cite Roussel2003.
 * @tparam DIM Dimension of the current simulation.
 * @tparam A Limiting convergence order of the overall scheme (space + time).
 */
template<Dimension DIM, unsigned int A>
class ThresholdComputer {

   static constexpr double D_     = static_cast<double>( DTI( DIM ) );//Capitalization according to paper
   static constexpr double alpha_ = static_cast<double>( A );
   unsigned int const maximum_level_;
   double const reference_epsilon_;

   /**
    * @brief Computes the reference epsilon used in the multiresolution analysis according to \cite Roussel2013.
    * @param user_epsilon_reference The threshold value the user wants to ensure on the given level.
    * @param user_level_of_reference The level at which the given threshold should be ensured.
    * @return Epsilon Reference to be used in multiresolution analysis.
    */
   inline double ComputeReferenceEpsilon( double const user_epsilon_reference, unsigned int const user_level_of_reference ) const {
      return user_epsilon_reference * std::pow( 2.0, -1.0 * ( alpha_ + 1.0 ) * ( double( maximum_level_ ) - double( user_level_of_reference ) ) );
   }

public:
   /**
    * @brief Standard constructor.
    * @param maximum_level Maximum level in the simulation.
    * @param user_reference_level Reference level provided by user input for the multiresolution threshold.
    * @param user_reference_epsilon Reference epsilon provided by user input for the multiresolution threshold.
    */
   explicit ThresholdComputer( unsigned int const maximum_level, unsigned int const user_reference_level, double const user_reference_epsilon ) : maximum_level_( maximum_level ),
                                                                                                                                                  reference_epsilon_( ComputeReferenceEpsilon( user_reference_epsilon, user_reference_level ) ) {
      // Empty besides initializer list
   }
   ThresholdComputer()                           = delete;
   ~ThresholdComputer()                          = default;
   ThresholdComputer( ThresholdComputer const& ) = default;
   ThresholdComputer& operator=( ThresholdComputer const& ) = delete;
   ThresholdComputer( ThresholdComputer&& )                 = default;
   ThresholdComputer& operator=( ThresholdComputer&& ) = delete;

   /**
    * @brief Gives the level-dependent threshold for the multiresolution analysis.
    * @param level The level on which the threshold is computed.
    * @return Computed threshold.
    */
   double ThresholdOnLevel( unsigned int const level ) const {
      return reference_epsilon_ * std::pow( 2.0, -1.0 * D_ * ( double( maximum_level_ ) - double( level ) ) );
   }
};

using Thresholder = ThresholdComputer<CC::DIM(), CC::STDO()>;

#endif// THRESHOLD_COMPUTER_H
