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
#include "input_output/utilities/xml_utilities.h"

#include <stdexcept>
#include "utilities/string_operations.h"

namespace XmlUtilities {

   namespace {
      /**
       * @brief Propagates errors tags from the current node, until the final stage.
       * @param node The node from which the propagated error is created.
       * @return The full error string.
       */
      std::string CreatePropagatedError( tinyxml2::XMLElement const* node ) {
         // Get the names of all parent nodes for this specific node
         std::string error_tags( "<" + std::string( node->Name() ) + ">" );
         // Get all parent nodes and add subsequently
         tinyxml2::XMLElement const* parent_node = node->Parent()->ToElement();
         while( parent_node != nullptr ) {
            std::string const parent_name( parent_node->Name() );
            error_tags += " under <" + parent_name + ">";
            parent_node = parent_node->Parent()->ToElement();
         }
         return error_tags;
      }

      /**
       * @brief Throws an invalid argument error with type from the current node, until the final stage.
       * @param node The node from which the propagated error is thrown.
       * @param type The type of the invalid argument that is thrown in the error message.
       */
      void ThrowPropagatedInvalidArgumentError( tinyxml2::XMLElement const* node, std::string const type ) {
         // Get the names of all parent nodes for this specific node
         std::string const error_tags( CreatePropagatedError( node ) );
         // Print the full error message
         throw std::invalid_argument( "Type error while reading argument ( " + type + " ) for " + error_tags + "!" );
      }

      /**
       * @brief Throws a logic error from the current node, until the final stage.
       * @param parent_node The node from which the propagated error is thrown.
       * @param child_name The additional child name that does not exist for the parent node.
       */
      void ThrowPropagatedLogicError( tinyxml2::XMLElement const* parent_node, std::string const child_name ) {
         // Get the names of all parent nodes for this specific node
         std::string const error_tags( CreatePropagatedError( parent_node ) );
         // Print the full error message
         throw std::logic_error( "Please specify block with tag <" + child_name + "> " + error_tags + "in Xml input file!" );
      }
   }// namespace

   /**
    * @brief Gives the child element for the given list for the given XML document.
    * @param xml_document the full xml document containing all xml information.
    * @param child_names List of child names (passed in subsequent).
    * @return Last child node in list (throw error if not existent).
    */
   tinyxml2::XMLElement const* GetChild( tinyxml2::XMLDocument const& xml_document, std::vector<std::string> child_names ) {
      if( !child_names.empty() ) {
         // Obtain the child name and erase it from the vector
         std::string const child_name( child_names.front() );
         child_names.erase( std::begin( child_names ) );

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
    * @brief Gives the child node for a given list of tags (iterate through list up to the end).
    * @param parent_node Parent node from which the child node list starts.
    * @param child_names List of child names (passed in subsequent).
    * @return Last child node in list (throw error if not existent).
    */
   tinyxml2::XMLElement const* GetChild( tinyxml2::XMLElement const* parent_node, std::vector<std::string> child_names ) {
      if( !child_names.empty() ) {
         // Obtain the child name and erase it from the vector
         std::string const child_name( child_names.front() );
         child_names.erase( std::begin( child_names ) );
         // Get the child node for the front element with the single string function (error back propagation)
         tinyxml2::XMLElement const* child_node = parent_node->FirstChildElement( child_name.c_str() );
         // if child does not exist throw error
         if( child_node == nullptr ) {
            ThrowPropagatedLogicError( parent_node, child_name );
         }
         // otherwise return the child of of the next child
         return GetChild( child_node, child_names );
      }
      // If this is the last name simply return it to finalize recursive call
      return parent_node;
   }

   /**
    * @brief Gives all child nodes of the given parent node that have the same given name.
    * @param parent_node The parent node.
    * @param child_name The child for which the nodes should be given.
    * @return The vector with all child nodes.
    */
   std::vector<tinyxml2::XMLElement const*> GetChilds( tinyxml2::XMLElement const* parent_node, std::string const& child_name ) {
      // The vector that is returned
      std::vector<tinyxml2::XMLElement const*> child_nodes;
      // Get the first child of the parent node
      tinyxml2::XMLElement const* child_node = parent_node->FirstChildElement();
      // Loop through all childs and store the nodes that correspond to the name
      while( child_node != nullptr ) {
         // Extract the name
         std::string const name = child_node->Name();
         if( child_node->Name() == child_name ) {
            child_nodes.push_back( child_node );
         }
         child_node = child_node->NextSiblingElement();
      }
      return child_nodes;
   }

   /**
    * @brief Checks whether the child exists for the given parent node.
    * @param xml_document the full xml document containing all xml information.
    * @param child_names List of child names (passed in subsequent).
    * @return True if exists.
    */
   bool ChildExists( tinyxml2::XMLDocument const& xml_document, std::vector<std::string> child_names ) {
      try {
         return GetChild( xml_document, child_names ) != nullptr;
      } catch( std::logic_error const& err ) {
         return false;
      }
   }

   /**
    * @brief Checks whether a child exists for the given parent node.
    * @param parent_node The parent node.
    * @param child_name The child_name.
    * @return true if child tag exists.
    */
   bool ChildExists( tinyxml2::XMLElement const* parent_node, std::string child_name ) {
      tinyxml2::XMLElement const* child_node = parent_node->FirstChildElement( child_name.c_str() );
      return child_node != nullptr;
   }

   /**
    * @brief Reads out a numeric value from an XML node, treats and converts it into a double value.
    * @param node The XML node holding the desired information.
    * @return The read-out and converted value.
    *
    * @note If the given type of the node does not coincide with a double value an error is thrown.
    */
   double ReadDouble( tinyxml2::XMLElement const* node ) {
      double value;
      if( node->QueryDoubleText( &value ) != tinyxml2::XML_SUCCESS ) {
         ThrowPropagatedInvalidArgumentError( node, "double" );
         throw 0;
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
   int ReadInt( tinyxml2::XMLElement const* node ) {
      int value;
      if( node->QueryIntText( &value ) != tinyxml2::XML_SUCCESS ) {
         ThrowPropagatedInvalidArgumentError( node, "int" );
         throw 0;
      } else {
         return value;
      }
   }

   /**
    * @brief Reads out a numeric value from an XML node, treats and converts it into an long int value.
    * @param node The XML node holding the desired information.
    * @return The read-out and converted value.
    *
    * @note If the given type of the node does not coincide with an integer value an error is thrown.
    */
   std::int64_t ReadLongInt( tinyxml2::XMLElement const* node ) {
      std::int64_t value;
      if( node->QueryInt64Text( &value ) != tinyxml2::XML_SUCCESS ) {
         ThrowPropagatedInvalidArgumentError( node, "long int" );
         throw 0;
      } else {
         return value;
      }
   }

   /**
    * @brief Reads out a numeric value from an XML node, treats and converts it into an unsigned int value.
    * @param node The XML node holding the desired information.
    * @return The read-out and converted value.
    *
    * @note If the given type of the node does not coincide with an unsigned integer value an error is thrown.
    */
   unsigned int ReadUnsignedInt( tinyxml2::XMLElement const* node ) {
      try {
         int const value = ReadInt( node );
         if( value < 0 ) {
            throw std::invalid_argument( "" );
         }
         return static_cast<unsigned int>( value );
      } catch( std::invalid_argument const& ) {
         ThrowPropagatedInvalidArgumentError( node, "unsigned int" );
         throw 0;
      }
   }

   /**
    * @brief Reads out a numeric value from an XML node, treats and converts it into an unsigned long int value.
    * @param node The XML node holding the desired information.
    * @return The read-out and converted value.
    *
    * @note If the given type of the node does not coincide with an unsigned integer value an error is thrown.
    */
   std::uint64_t ReadUnsignedLongInt( tinyxml2::XMLElement const* node ) {
      try {
         std::int64_t const value = ReadLongInt( node );
         if( value < 0 ) {
            throw std::invalid_argument( "" );
         }
         return static_cast<unsigned int>( value );
      } catch( std::invalid_argument const& ) {
         ThrowPropagatedInvalidArgumentError( node, "unsigned long int" );
         throw 0;
      }
   }

   /**
    * @brief Reads out a string value from an XML node.
    * @param node The XML node holding the desired information.
    * @return Found String
    *
    * @note If the given type of the node does not coincide with a string an error is thrown.
    */
   std::string ReadString( tinyxml2::XMLElement const* node ) {
      // Read the text from the node
      char const* text = node->GetText();
      std::string const value( text != nullptr ? StringOperations::Trim( text ) : "" );
      if( value.empty() || !value.compare( "" ) ) {
         ThrowPropagatedInvalidArgumentError( node, "string" );
         throw 0;
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
   std::vector<double> ReadTimeStamps( tinyxml2::XMLElement const* parent_node ) {
      // Base name for all tags
      std::string timestamp_name;
      std::string base_name( "ts" );
      // vector where the final time stamps are written into
      std::vector<double> timestamps;
      // Read until final number is reached
      bool keep_reading = true;
      int index         = 1;
      while( keep_reading ) {
         timestamp_name                   = base_name + std::to_string( index++ );
         tinyxml2::XMLElement const* node = parent_node->FirstChildElement( timestamp_name.c_str() );
         if( node == nullptr ) {
            keep_reading = false;
         } else {
            timestamps.emplace_back( ReadDouble( node ) );
         }
      }

      return timestamps;
   }
}// namespace XmlUtilities
