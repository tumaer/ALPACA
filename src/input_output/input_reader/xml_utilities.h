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
#ifndef XML_UTILITIES_H
#define XML_UTILITIES_H

#include <stdexcept>
#include <vector>
#include <iostream>

#include <tinyxml2.h>

/**
 * @brief In this name space all functions are defiend that simplify the reading process of xml input files.
 */
namespace XmlUtilities {

   /**
    * @brief Gives the child node for a given list of tags (iterate through list up to the end).
    * @param parent_node Parent node from which the child node list starts.
    * @param child_names List of child names (passed in subsequent).
    * @return Last child node in list (throw error if not existent).
    */
   inline tinyxml2::XMLElement const* GetChild( tinyxml2::XMLElement const *parent_node, std::vector<std::string> child_names  ) {
      if( !child_names.empty() ) {
         // Obtain the child name and erase it from the vector
         std::string const child_name( child_names.front() );
         child_names.erase( child_names.begin() );
         // try to get the child node and proceed with remaining childs
         try {
            // Get the child node for the front element with the single string function (error back propagation)
            tinyxml2::XMLElement const* child_node = parent_node->FirstChildElement( child_name.c_str() );
            // if child does not exist throw error
            if( child_node == nullptr ) {
               throw std::logic_error( "Please specify block with tag <" + child_name + ">" );
            }
            // otherwise return the child of of the next child
            return GetChild( child_node, child_names );

         } catch ( std::logic_error const& err ) {
            std::string const parent_name( parent_node->Name() );
            // Back propagate the error message to its original parent node
            throw std::logic_error( std::string( err.what() ) + " under <" + parent_name + ">" );
         }
      }
      // If this is the last name simply return it to finalize recursive call
      return parent_node;
   }

   /**
    * @brief Gives the child element for the given list for the given XML document.
    * @param xml_document the full xml document containing all xml information.
    * @param child_names List of child names (passed in subsequent).
    * @return Last child node in list (throw error if not existent).
    */
   inline tinyxml2::XMLElement const* GetChild( tinyxml2::XMLDocument const& xml_document, std::vector<std::string> child_names  ) {
      if( !child_names.empty() ) {
         // Obtain the child name and erase it from the vector
         std::string const child_name( child_names.front() );
         child_names.erase( child_names.begin() );

         // Get the child node for the front element with the single string function (error back propagation)
         tinyxml2::XMLElement const* child_node = xml_document.FirstChildElement( child_name.c_str() );
         // if child does not exist throw error
         if( child_node == nullptr ) {
            throw std::logic_error( "Please specify block with tag <" + child_name + "> in Xml input file!" );
         }
         // otherwise return the child of of the next child (call XMl element function)
         return GetChild( child_node, child_names );

      } else {
         throw std::logic_error( "You must specify at least one child that should be obtained from the Xml-Document!" );
      }
   }

   /**
    * @brief Reads out a numeric value from an XML node, treats and converts it into a double value.
    * @param node The XML node holding the desired information.
    * @return The read-out and converted value.
    *
    * @note If the given type of the node does not coincide with a double value an error is thrown.
    */
   inline double ReadDouble( tinyxml2::XMLElement const *node ) {
      double value;
      if( node->QueryDoubleText( &value ) != tinyxml2::XML_SUCCESS ) {
         // Get the names of all parent nodes for this specific node
         std::string error_tags( "<" + std::string( node->Name() ) + ">" );
         // Get all parent nodes and add subsequently
         tinyxml2::XMLElement const* parent_node = node->Parent()->ToElement();
         while( parent_node != nullptr ) {
            std::string const parent_name( parent_node->Name() );
            error_tags += " under <" + parent_name + ">";
            parent_node = parent_node->Parent()->ToElement();
         }
         // Print the full error message
         throw std::invalid_argument( "Type error while reading argument ( double ) for " + error_tags + "!" );
      } else {
         return value;
      }
   }

   /**
    * @brief Reads out a numeric value from an XML node, treats and converts it into an int value.
    * @param node The XML node holding the desired information.
    * @return The read-out and converted value.
    *
    * @note If the given type of the node does not coincide with an integer value an error is thrown.
    */
   inline int ReadInt( tinyxml2::XMLElement const *node ) {
      int value;
      if( node->QueryIntText( &value ) != tinyxml2::XML_SUCCESS ) {
         // Get the names of all parent nodes for this specific node
         std::string error_tags( "<" + std::string( node->Name() ) + ">" );
         // Get all parent nodes and add subsequently
         tinyxml2::XMLElement const* parent_node = node->Parent()->ToElement();
         while( parent_node != nullptr ) {
            std::string const parent_name( parent_node->Name() );
            error_tags += " under <" + parent_name + ">";
            parent_node = parent_node->Parent()->ToElement();
         }
         // Print the full error message
         throw std::invalid_argument( "Type error while reading argument ( int ) for " + error_tags + "!" );
      } else {
         return value;
      }
   }

   /**
    * @brief Reads out a string value from an XML node.
    * @param node The XML node holding the desired information.
    * @return Found String
    *
    * @note If the given type of the node does not coincide with a string an error is thrown.
    */
   inline std::string ReadString( tinyxml2::XMLElement const *node ) {

      std::string const value( node->GetText() );

      if( value.empty() || !value.compare( "" ) ) {
         // Get the names of all parent nodes for this specific node
         std::string error_tags( "<" + std::string( node->Name() ) + ">" );
         // Get all parent nodes and add subsequently
         tinyxml2::XMLElement const* parent_node = node->Parent()->ToElement();
         while( parent_node != nullptr ) {
            std::string const parent_name( parent_node->Name() );
            error_tags += " under <" + parent_name + ">";
            parent_node = parent_node->Parent()->ToElement();
         }
         // Print the full error message
         throw std::invalid_argument( "Type error while reading argument ( string ) for " + error_tags + "!" );
      }
      return value;
   }

   /**
    * @brief Reads out all time stamps of a given parent node.
    * @param parent_node Where the time stamps should be read from.
    * @return All read time stamps.
    *
    * @note If the order starting with 1 is broken, following time steps are not read anymore.
    */
   inline std::vector<double> ReadTimeStamps( tinyxml2::XMLElement const *parent_node ) {
      // Base name for all tags
      std::string timestamp_name;
      std::string base_name( "ts" );
      // vector where the final time stamps are written into
      std::vector<double> timestamps;
      // Read until final number is reached
      bool keep_reading = true;
      int index = 1;
      while( keep_reading ) {
         timestamp_name = base_name + std::to_string( index++ );
         tinyxml2::XMLElement const* node = parent_node->FirstChildElement( timestamp_name.c_str() );
         if( node == nullptr ) {
            keep_reading = false;
         } else {
            timestamps.emplace_back( ReadDouble( node ) );
         }
      }

      return timestamps;
   }
}

#endif // XML_UTILITIES_H