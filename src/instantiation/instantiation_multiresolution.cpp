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
#include "instantiation/instantiation_multiresolution.h"

namespace Instantiation {

   /**
    * @brief Creates a Multiresolution object.
    * @param input_reader Reader that provides access to the full data of the input file.
    * @param topology_manager Class providing global (on all ranks) node information.
    * @return A multiresolution object with thresholding conditions accroding to given input.
    */
   Multiresolution InstantiateMultiresolution( InputReader const& input_reader, TopologyManager const& topology_manager ) {

      // Get the data local (for logging)
      unsigned int const maximum_level           = topology_manager.GetMaximumLevel();
      unsigned int const epsilon_reference_level = input_reader.GetMultiResolutionReader().ReadEpsilonLevelReference();
      double const epsilon_reference             = input_reader.GetMultiResolutionReader().ReadEpsilonReference();

      // Create the thresholder
      Thresholder thresholder = Thresholder( maximum_level, epsilon_reference_level, epsilon_reference );

      // Log data
      LogWriter& logger = LogWriter::Instance();
      logger.LogMessage( " " );
      if( maximum_level > 1 ) {
         // tmp string for maximum size determination
         std::string const level_string( std::to_string( maximum_level - 1 ) + std::to_string( maximum_level ) );

         logger.LogMessage( "Epsilon Level 0 and Level 1" + std::string( level_string.size() - 2, ' ' ) + " : " + StringOperations::ToScientificNotationString( thresholder.ThresholdOnLevel( 1 ), 9 ) );
         logger.LogMessage( "Epsilon Level " + std::to_string( maximum_level - 1 ) + " and Level " + std::to_string( maximum_level ) + " : " + StringOperations::ToScientificNotationString( thresholder.ThresholdOnLevel( maximum_level ), 9 ) );
      } else if( maximum_level == 1 ) {
         logger.LogMessage( "Homogenous Mesh - Level 1 cannot be coarsened. Provided Epsilon is ignored" );
      } else {
         logger.LogMessage( "Homogenous Mesh - Provided Epsilon is ignored " );
      }
      logger.LogMessage( " " );

      // Return the created multiresolution with the thresholder
      return Multiresolution( std::move( thresholder ) );
   }
}// namespace Instantiation
