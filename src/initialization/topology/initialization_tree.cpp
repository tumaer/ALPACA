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
#include "initialization/topology/initialization_tree.h"

namespace Initialization {

   /**
    * @brief Initializes the complete tree class with the given input classes
    * @param input_reader Reader that provides access to the full data of the input file
    * @param topology_manager Class providing global (on all ranks) node information  
    * @param unit_handler Instance to provide (non-)dimensionalization of values
    * @return The fully initialized tree class 
    */
   Tree InitializeTree( InputReader const& input_reader, TopologyManager & topology_manager, UnitHandler const& unit_handler ) {

      // Get the data for initializing and/or logging 
      unsigned int const maximum_level = topology_manager.GetMaximumLevel();
      std::array<unsigned int, 3> const number_of_nodes = topology_manager.GetNumberOfNodesOnLevelZero();
      double const node_size_on_level_zero = input_reader.GetMultiResolutionReader().ReadNodeSizeOnLevelZero();

      // Logging of data 
      LogWriter & logger = LogWriter::Instance();
      std::string tmp_string; 

      logger.LogMessage( " " );
      // Domain size 
      tmp_string =      "Domain size                : "   + StringOperations::ToScientificNotationString( number_of_nodes[0] * node_size_on_level_zero, 6 );
      tmp_string += CC::DIM() != Dimension::One ?   " x " + StringOperations::ToScientificNotationString( number_of_nodes[1] * node_size_on_level_zero, 6 ) : "";
      tmp_string += CC::DIM() == Dimension::Three ? " x " + StringOperations::ToScientificNotationString( number_of_nodes[2] * node_size_on_level_zero, 6 ) : "";
      logger.LogMessage( tmp_string );
      // Internal cells per block 
      logger.LogMessage( "Internal cell per block    : "  + std::to_string( CC::ICX() ) );
      // Maximum level 
      logger.LogMessage( "Maximum level              : "  + std::to_string( maximum_level ) );
      // Resolution level zero 
      tmp_string =       "Resolution on level zero   : "  + std::to_string( number_of_nodes[0] * CC::ICX() );
      tmp_string += CC::DIM() != Dimension::One ?   " x " + std::to_string( number_of_nodes[1] * CC::ICX() ) : "";
      tmp_string += CC::DIM() == Dimension::Three ? " x " + std::to_string( number_of_nodes[2] * CC::ICX() ) : "";
      logger.LogMessage( tmp_string + " internal cells" );
      // Cell size on level zero 
      logger.LogMessage( "Cell size on level zero    : "        
                        + StringOperations::ToScientificNotationString( node_size_on_level_zero / double( CC::ICX() ), 9 ) );
      // Resolution level maximum 
      tmp_string =       "Resolution on maximum level: "  + std::to_string( number_of_nodes[0] * CC::ICX() * ( 1 << maximum_level ) );
      tmp_string += CC::DIM() != Dimension::One ?   " x " + std::to_string( number_of_nodes[1] * CC::ICX() * ( 1 << maximum_level ) ) : "";
      tmp_string += CC::DIM() == Dimension::Three ? " x " + std::to_string( number_of_nodes[2] * CC::ICX() * ( 1 << maximum_level ) ) : "";
      logger.LogMessage( tmp_string + " internal cells" );
      // Cell size on level maximum 
      logger.LogMessage( "Cell size on maximum level : " 
                        + StringOperations::ToScientificNotationString( node_size_on_level_zero / double( CC::ICX() ) / double( 1 << maximum_level ), 9 ) );
      logger.LogMessage( " " );

      // return the initialized tree 
      return Tree( topology_manager, 
                   topology_manager.GetMaximumLevel(), 
                   unit_handler.NonDimensionalizeValue( node_size_on_level_zero, UnitType::Length ) );

   }
}