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
#include "instantiation/input_output/instantiation_input_reader.h"

#include <stdexcept>
#include <memory>

#include <tinyxml2.h>
#include "input_output/input_reader/input_definitions.h"
#include "utilities/file_operations.h"

#include "input_output/input_reader/boundary_condition_reader/xml_boundary_condition_reader.h"
#include "input_output/input_reader/initial_condition_reader/xml_initial_condition_reader.h"
#include "input_output/input_reader/material_reader/xml_material_reader.h"
#include "input_output/input_reader/multi_resolution_reader/xml_multi_resolution_reader.h"
#include "input_output/input_reader/dimensionalization_reader/xml_dimensionalization_reader.h"
#include "input_output/input_reader/output_reader/xml_output_reader.h"
#include "input_output/input_reader/restart_reader/xml_restart_reader.h"
#include "input_output/input_reader/source_term_reader/xml_source_term_reader.h"
#include "input_output/input_reader/time_control_reader/xml_time_control_reader.h"

namespace Instantiation {

   /**
    * @brief Instantiates the full input reader class with the given input file.
    * @param input_filename Name of the file use for input.
    * @return The fully instantiated InputReader class.
    */
   InputReader InstantiateInputReader( std::string const& input_filename ) {
      // Check whether the file exists
      if( !FileOperations::CheckIfPathExists( input_filename ) ) {
         throw std::logic_error( "Input file " + input_filename + " does not exist!" );
      }

      // Determine the input type
      InputType const input_type( StringToInputType( FileOperations::GetFileExtension( input_filename ) ) );
      // Instantiate correct reader
      switch( input_type ) {
         case InputType::Xml: {
            // Open the file (here std::make_shared not possible)
            // shared pinter required to distribute the open input file on different reader
            std::shared_ptr<tinyxml2::XMLDocument> input_file( new tinyxml2::XMLDocument );
            tinyxml2::XMLError error = input_file->LoadFile( input_filename.c_str() );
            // Check if eversthing worked properly
            if( error != tinyxml2::XML_SUCCESS ) {
               throw std::logic_error( "Syntax error parsing the XML inputfile file, check opening and closing tags!" );
            }
            // Create the input reader properly
            return InputReader( input_filename, input_type,
                                std::make_unique<XmlMaterialReader const>( input_file ),
                                std::make_unique<XmlBoundaryConditionReader const>( input_file ),
                                std::make_unique<XmlInitialConditionReader const>( input_file ),
                                std::make_unique<XmlMultiResolutionReader const>( input_file ),
                                std::make_unique<XmlDimensionalizationReader const>( input_file ),
                                std::make_unique<XmlOutputReader const>( input_file ),
                                std::make_unique<XmlRestartReader const>( input_file ),
                                std::make_unique<XmlSourceTermReader const>( input_file ),
                                std::make_unique<XmlTimeControlReader const>( input_file ) );
         }

         default: {
            throw std::invalid_argument( "Input file extension not known! Cannot choose correct reader!" );
         }
      }
   }
}// namespace Instantiation
