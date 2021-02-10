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

#include <string>
#include <memory>

#include "input_output/input_reader/dimensionalization_reader/dimensionalization_reader.h"
#include "input_output/input_reader/dimensionalization_reader/xml_dimensionalization_reader.h"

SCENARIO( "Check that the xml dimensionalization reader works properly", "[1rank]" ) {
   GIVEN( "A xml document with the valid content to read the dimensionalization data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <dimensionalization>"
                                  "     <lengthReference> 1.0 </lengthReference>"
                                  "     <velocityReference> 2.0 </velocityReference>"
                                  "     <densityReference>   3.0 </densityReference>"
                                  "     <temperatureReference>   4.0 </temperatureReference>"
                                  "  </dimensionalization>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<DimensionalizationReader const> const reader( std::make_unique<XmlDimensionalizationReader const>( xml_tree ) );
      WHEN( "The reference length is read." ) {
         THEN( "The length should be 1.0" ) {
            REQUIRE( reader->ReadReferenceLength() == 1.0 );
         }
      }
      WHEN( "The reference velocity is read." ) {
         THEN( "The velocity should be 2.0" ) {
            REQUIRE( reader->ReadReferenceVelocity() == 2.0 );
         }
      }
      WHEN( "The reference density is read." ) {
         THEN( "The density should be 3.0" ) {
            REQUIRE( reader->ReadReferenceDensity() == 3.0 );
         }
      }
      WHEN( "The reference temperature is read." ) {
         THEN( "The temperature should be 4.0" ) {
            REQUIRE( reader->ReadReferenceTemperature() == 4.0 );
         }
      }
   }
   GIVEN( "A xml document with invalid content to read the dimensionalization data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <dimensionalization>"
                                  "     <lengthReference> -1.0 </lengthReference>"
                                  "     <velocityReference> -2.0 </velocityReference>"
                                  "     <densityReference>  - 3.0 </densityReference>"
                                  "     <temperatureReference>   -4.0 </temperatureReference>"
                                  "  </dimensionalization>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<DimensionalizationReader const> const reader( std::make_unique<XmlDimensionalizationReader const>( xml_tree ) );
      WHEN( "Data is read from the ree." ) {
         THEN( "An std::invalid_argument exception should be thrown" ) {
            REQUIRE_THROWS_AS( reader->ReadReferenceLength(), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadReferenceVelocity(), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadReferenceDensity(), std::invalid_argument );
            REQUIRE_THROWS_AS( reader->ReadReferenceTemperature(), std::invalid_argument );
         }
      }
   }

   GIVEN( "A xml document with the non-existent tags to read the dimensionalization data." ) {
      std::string const xml_data( "<configuration>"
                                  "  <dimensionalization>"
                                  "  </dimensionalization>"
                                  "</configuration>" );
      // Create the xml document
      std::shared_ptr<tinyxml2::XMLDocument> xml_tree( std::make_shared<tinyxml2::XMLDocument>() );
      xml_tree->Parse( xml_data.c_str() );
      // Create the xml reader
      std::unique_ptr<DimensionalizationReader const> const reader( std::make_unique<XmlDimensionalizationReader const>( xml_tree ) );
      WHEN( "Any data is read." ) {
         THEN( "All should throw std::logic_error exception" ) {
            REQUIRE_THROWS_AS( reader->ReadReferenceLength(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadReferenceVelocity(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadReferenceDensity(), std::logic_error );
            REQUIRE_THROWS_AS( reader->ReadReferenceTemperature(), std::logic_error );
         }
      }
   }
}